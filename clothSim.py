from __future__ import division
import time
import math
import numpy as np
import maya.cmds as cmds
from collections import namedtuple

class ParticleSystem:
    def __init__(self, sub_x=6, sub_y=6, len_x=6, len_y=6, pivot_loc=(0,0,0), num_iter=10, t_step=0.05):

        self.m_time_step_num = 1
        self.iter_count = 0

        self.m_row_vtx = sub_x+1
        self.m_col_vtx = sub_y+1

        self.m_cur_vtx_pos = []
        self.m_vtx_clus = []
        self.__create_cloth_mesh(sub_x, sub_y, len_x, len_y, pivot_loc) # name = m_cloth_mesh

        self.m_old_vtx_pos = np.copy(self.m_cur_vtx_pos)
        self.m_force_acc = [np.zeros(3)]*len(self.m_cur_vtx_pos)    # Force Accumulators
        self.m_constraints = []
        self.__create_all_constraints()

        self.fixed_vtx_pos_01 = np.copy(self.m_cur_vtx_pos[0])
        self.fixed_vtx_pos_02 = np.copy(self.m_cur_vtx_pos[sub_x])
        self.m_num_con_iterations = num_iter
        self.m_gravity = np.array([0.0,-9.8,0.0])
        self.m_time_step = np.float(t_step)

        self.m_sphere_colliders = {}


    def __create_cloth_mesh(self, sub_x, sub_y, len_x, len_y, pivot_loc):
        cmds.polyPlane(name='m_cloth_mesh', sx=sub_x, sy=sub_y, w=len_x, h=len_y)
 
        clus_grp = cmds.group(em=True, n='vtx_cluster_GRP')
        mesh_vtx = cmds.ls('m_cloth_mesh.vtx[:]', fl=True)
        for vtx in mesh_vtx:
            cur_clus = cmds.cluster(vtx)
            cmds.setAttr('%sHandle.visibility' % cur_clus[0], False)
            cmds.parent(cur_clus, clus_grp)
            cmds.setAttr('%sHandle.translate' % cur_clus[0], pivot_loc[0], pivot_loc[1], pivot_loc[2])
            self.m_vtx_clus.append(cur_clus[1])
        self.__get_cloth_vtx_pos()


    def __get_cloth_vtx_pos(self):
        cloth_vtx = cmds.xform('m_cloth_mesh.vtx[*]', q=True, ws=True, t=True)
        i=0
        while i<len(cloth_vtx):
            self.m_cur_vtx_pos.append(np.array(cloth_vtx[i:i+3]))
            i+=3

    def __create_all_constraints(self):
        num_vtx = (self.m_row_vtx)*(self.m_col_vtx)
        for i in range(num_vtx):
            self.__create_vtx_constraints(i)


    def __create_vtx_constraints(self, i):
        num_vtx = (self.m_row_vtx)*(self.m_col_vtx)
        # Right
        temp_vtx = i+1
        if num_vtx>temp_vtx:    # if a vertex to the right of 'i' exists
            if math.ceil((i+1)/self.m_row_vtx) == math.ceil((temp_vtx+1)/self.m_row_vtx):    # if 'i' and 'i+1' are in the same row
                self.m_constraints.append(Constraint(i,temp_vtx,self.__get_rest_length(i,temp_vtx)))
        # Rightx2
        temp_vtx = i+2
        if num_vtx>temp_vtx:
            if math.ceil((i+1)/self.m_row_vtx) == math.ceil((temp_vtx+1)/self.m_row_vtx):
                self.m_constraints.append(Constraint(i,temp_vtx,self.__get_rest_length(i,temp_vtx)))
        # Down and Diagonal
        temp_vtx = i+self.m_row_vtx
        if num_vtx>temp_vtx:    # if a vertex directly under 'i' exists
            self.m_constraints.append(Constraint(i,temp_vtx,self.__get_rest_length(i,temp_vtx)))

            diag_vtx = temp_vtx+1
            if num_vtx>diag_vtx:
                if math.ceil((temp_vtx+1)/self.m_row_vtx) == math.ceil((diag_vtx+1)/self.m_row_vtx):
                    self.m_constraints.append(Constraint(i,diag_vtx,self.__get_rest_length(i,diag_vtx)))

            diag_vtx = temp_vtx-1
            if math.ceil((temp_vtx+1)/self.m_row_vtx) == math.ceil((diag_vtx+1)/self.m_row_vtx):
                self.m_constraints.append(Constraint(i,diag_vtx,self.__get_rest_length(i,diag_vtx)))
        # Downx2
        temp_vtx = i+(self.m_row_vtx*2)
        if num_vtx>temp_vtx:    # if a vertex directly under 'i' exists
            self.m_constraints.append(Constraint(i,temp_vtx,self.__get_rest_length(i,temp_vtx)))


    def __get_rest_length(self, vtx_01, vtx_02):    # vtx_01 and vtx_02 are indices in 'ParticleSystem.m_cur_vtx_pos' with positions for these vertices.
        return np.linalg.norm(self.m_cur_vtx_pos[vtx_02]-self.m_cur_vtx_pos[vtx_01])


    def __accumulate_forces(self):
        for i in range(len(self.m_cur_vtx_pos)):
            self.m_force_acc[i] = np.copy(self.m_gravity)


    def __perform_verlet(self):    # Verlet Integration Step
        for i in range(len(self.m_cur_vtx_pos)):
            x = np.copy(self.m_cur_vtx_pos[i])
            temp = np.copy(x)
            old_x = np.copy(self.m_old_vtx_pos[i])
            a = np.copy(self.m_force_acc[i])
            x += x-old_x+a*self.m_time_step*self.m_time_step
            old_x = temp
            self.m_cur_vtx_pos[i] = np.copy(x)
            self.m_old_vtx_pos[i] = np.copy(old_x)


    def __satisfy_constraints(self):
        for k in range(self.m_num_con_iterations):
            
            self.__satisfy_sphere_collider_constraints()
            # self.__satisfy_box_collider_constraints()
            
            for i in range(len(self.m_constraints)):
                c = self.m_constraints[i]

                x1 = self.m_cur_vtx_pos[c.particle_a]
                x2 = self.m_cur_vtx_pos[c.particle_b]
                delta = x2-x1
                deltalength = math.sqrt(np.dot(delta, delta))
                diff = (deltalength-c.rest_length)/deltalength * 0.1
                x1 += delta*0.5*diff
                x2 -= delta*0.5*diff
                # self.m_cur_vtx_pos[c.particle_a] = x1
                # self.m_cur_vtx_pos[c.particle_b] = x2  

            # if self.iter_count<90:
            self.m_cur_vtx_pos[0] = np.copy(self.fixed_vtx_pos_01)
            self.m_cur_vtx_pos[self.m_row_vtx-1] = np.copy(self.fixed_vtx_pos_02)

        # self.iter_count+=1

    def time_step(self):
        self.__accumulate_forces()
        self.__perform_verlet()
        self.__satisfy_constraints()
        self.__set_vtx_pos()
        self.__set_keyframe()


    def __set_vtx_pos_old(self):
        num_vtx = (self.m_row_vtx)*(self.m_col_vtx)
        for i in range(num_vtx):
            cmds.xform('m_cloth_mesh.vtx[%d]' % i, ws=True, t=self.m_cur_vtx_pos[i])


    def __set_vtx_pos(self):
        num_vtx = (self.m_row_vtx)*(self.m_col_vtx)
        for i in range(num_vtx):
            clus_offset = cmds.getAttr('%s.origin' % self.m_vtx_clus[i])
            delta = (self.m_cur_vtx_pos[i]-clus_offset)[0]            
            cmds.xform('%s' % self.m_vtx_clus[i], ws=True, t=delta)

    def __set_keyframe(self):
        cmds.currentTime(self.m_time_step_num)
        cmds.select(cmds.listRelatives('vtx_cluster_GRP'))
        cmds.setKeyframe(attribute='translate')
        self.m_time_step_num+=1

    def add_sphere_collider(self, rad=1):
        cmds.polySphere(radius=rad)
        col_name = cmds.ls(sl=True)[0]    # col = Collider
        self.m_sphere_colliders[col_name]= rad


    def del_sphere_collider(self):
        del_list = cmds.ls(sl=True)
        for item in del_list:
            cmds.delete(item)
            del self.m_sphere_colliders[item]


    def __satisfy_sphere_collider_constraints(self):
        for cur_s_col in self.m_sphere_colliders:
            cur_s_pivot = cmds.getAttr('%s.translate' % cur_s_col)[0]
            cur_s_r_sqr = math.pow(self.m_sphere_colliders[cur_s_col],2)
            i=0
            for cur_particle in self.m_cur_vtx_pos:
                x = math.pow((cur_particle[0]-cur_s_pivot[0]),2)
                y = math.pow((cur_particle[1]-cur_s_pivot[1]),2)
                z = math.pow((cur_particle[2]-cur_s_pivot[2]),2)
                if x+y+z<cur_s_r_sqr:    # if 'cur_particle' is inside sphere
                    self.__handle_collisions(cur_particle, cur_s_pivot, self.m_sphere_colliders[cur_s_col], i)
                i+=1


    def __handle_collisions(self, x, p, r, i):
        self.m_cur_vtx_pos[i] = (((x-p)/np.linalg.norm(x-p))*(1.025*r))+p

    
    def __satisfy_box_collider_constraints(self):
        for i in range(len(self.m_cur_vtx_pos)):
            x = self.m_cur_vtx_pos[i]
            x = np.minimum(np.maximum(x, np.array([-10.0,0.1,-10.0])), np.array([100.0,100.0,100.0]))
            self.m_cur_vtx_pos[i] = x
        # print x
        

class Constraint:
    def __init__(self, particle_a, particle_b, rest_length):
        self.particle_a = particle_a    # 'ParticleSystem.m_cloth_vtx' index with the particle
        self.particle_b = particle_b
        self.rest_length = rest_length  


cloth_sim = ParticleSystem(sub_x=24, sub_y=24, len_x=8, len_y=8, pivot_loc=(0,20,0), num_iter=50, t_step=0.0125)

cloth_sim.add_sphere_collider(2)
# cloth_sim.del_sphere_collider()
# cloth_sim.time_step()

'''
t = time.time()
for i in range(350):
    cloth_sim.time_step()
print 'time: %s' %(time.time() - t)

sim_name = "cons_verlet"
sim_num = 1
cmds.playblast(f="%s_%d.mv" % (sim_name, sim_num), v=False)
sim_num+=1
'''