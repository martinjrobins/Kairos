/* 
 * nsvc.cpp
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of Smoldyn.
 *
 * Smoldyn is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Smoldyn is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Smoldyn.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Apr 16, 2013
 *      Author: mrobins
 */

#include "nsvc.h"
#include "smoldynfuncs.h"
#include <sstream>
#include <set>

void nsv_init() {
	Kairos::init(0,NULL);
}

NextSubvolumeMethod* nsv_new(double* min, double* max, double* dx, int n) {
	using namespace Kairos;
	Vect3d min_vect(0.0,0.0,0.0);
	Vect3d max_vect(1.0,1.0,1.0);
	Vect3d dx_vect(1.0,1.0,1.0);
	for (int i = 0; i < n; ++i) {
		min_vect[i] = min[i];
		max_vect[i] = max[i];
		dx_vect[i] = dx[i];
	}
	StructuredGrid *grid = new StructuredGrid(min_vect,max_vect,dx_vect);
	NextSubvolumeMethod* nsv = new NextSubvolumeMethod(*grid);
	return nsv;
}

void nsv_delete(NextSubvolumeMethod* nsv) {
	delete &nsv->get_grid();
	delete nsv;
}

void nsv_print(NextSubvolumeMethod* nsv, char* buffer) {
	std::ostringstream ss;
	ss << std::endl << *nsv << std::endl;
	ss.str().copy(buffer,ss.str().length(),0);
	buffer[ss.str().length()] = '\0';
}

void nsv_add_interface(NextSubvolumeMethod* nsv,double dt, double *start,double *end,double *norm,int dim) {
	using namespace Kairos;
	Vect3d min(0.0,0.0,0.0);
	Vect3d max(1.0,1.0,1.0);
	for (int i = 0; i < dim; ++i) {
		min[i] = start[i];
		max[i] = end[i];
	}
	if ((norm[0] == 1) || (norm[0] == -1)) {
		xrect interface(min,max,norm[0]);
		nsv->set_interface(interface,dt);

	} else if ((norm[1] == 1) || (norm[1] == -1)) {
		yrect interface(min,max,norm[1]);
		nsv->set_interface(interface,dt);

	} else if ((norm[2] == 1) || (norm[2] == -1)) {
		zrect interface(min,max,norm[2]);
		nsv->set_interface(interface,dt);

	}
}
void nsv_add_species(NextSubvolumeMethod* nsv,int id,double D) {
	using namespace Kairos;

	//TODO: assumes that species has not been added before
	Species *ns = new Species(id,D,&nsv->get_grid());
	nsv->add_diffusion(*ns);

}
void nsv_add_reaction(NextSubvolumeMethod* nsv,double rate,
						int nreactants,int *reactant_ids,
						int nproducts,int *product_ids) {

	//TODO: assumes that species has already been added via nsv_add_species
	using namespace Kairos;


	ReactionSide lhs;
	for (int i = 0; i < nreactants; ++i) {
		lhs = lhs + *(nsv->get_species(reactant_ids[i]));
	}
	ReactionSide rhs;
	for (int i = 0; i < nproducts; ++i) {
		rhs = rhs + *(nsv->get_species(product_ids[i]));
	}
	nsv->add_reaction(rate, lhs >> rhs);
}

void nsv_integrate(NextSubvolumeMethod* nsv,double dt, portstruct *port) {
	using namespace Kairos;
	//std::cout << "running lattice dt"<<std::endl;
	nsv->integrate(dt);
	//std::cout << " lattice t = "<<nsv->get_time()<<" smoldyn time = "<<port->portss->sim->time<<std::endl;

	/*
	 * put new particles in port
	 */
	simptr sim = port->portss->sim;

	const int ns = nsv->get_diffusing_species().size();
	int n = 0;
	for (int i = 0; i < ns; ++i) {
		n += nsv->get_diffusing_species()[i]->particles.size();
	}
	int *species = new int[n];
	double **positions = new double*[n];

	n = 0;
	for (int i = 0; i < ns; ++i) {
		Species *s = nsv->get_diffusing_species()[i];
		const int np = s->particles.size();
		for (int j = 0; j < np; ++j) {
			species[n] = s->id;
			positions[n] = s->particles[j].data();
			n++;
		}
		s->particles.clear();
	}
	int er = portputmols(sim,port,n,species[0],species,positions);
	if (er != 0) simLog(NULL,11,"ERROR: failure in nsv_integrate\n");

	delete species;
	delete positions;


	/*
	 * look in port for new particles
	 */

	//n = portgetmols2(sim,port,-1,MSall,0,positions);
	n = sim->mols->nl[port->llport];
	species = new int[n];
	positions = new double*[n];
	//n = portgetmols2(sim,port,-1,MSall,1,positions);

	std::set<int> dirty_indicies;
	int nout = 0;
	for (int i = 0; i < n; ++i) {
		moleculeptr m = sim->mols->live[port->llport][i];

		//check if already gone back through port
		double *crsspt = new double[sim->dim];
		enum PanelFace face1,face2;
		bool in = true;
		for (int j = 0; j < port->srf->npanel[PSrect]; ++j) {
			if (lineXpanel(m->posx,m->pos,port->srf->panels[PSrect][j],sim->dim,crsspt,&face1,&face2,NULL,NULL,NULL)) {
				if (face2==PFfront) {
					in = false;
				}
			}
//			if (ptinpanel(m->pos,port->srf->panels[j],sim->dim)) {
//				in_port = 1;
//				break;
//			}
		}
		delete crsspt;

		if (!in) {
			species[nout] = m->ident;
			positions[nout] = m->pos;
			nout++;
			continue;
		}

		Vect3d newr(0.5,0.5,0.5);
		for (int d = 0; d < sim->dim; ++d) {
			//assume reflective bc if outside domain (for now)
			double low = nsv->get_grid().get_low()[d];
			double high = nsv->get_grid().get_high()[d];
			if (m->pos[d] < low) {
				newr[d] = 2*low - m->pos[d];
			} else if (m->pos[d] > high) {
				newr[d] = 2*high - m->pos[d];
			} else {
				newr[d] = m->pos[d];
			}
		}

		const int ci = nsv->get_grid().get_cell_index(newr);
		nsv->get_species(m->ident)->copy_numbers[ci]++;
		dirty_indicies.insert(ci);
	}
	for (auto i : dirty_indicies) {
		nsv->recalc_priority(i);
	}


	//delete particles in port
	er = portgetmols(sim,port,-1,MSall,1);
	if (er != n) simLog(NULL,11,"ERROR: failure in nsv_integrate (while deleting particles from port)\n");

	//put particles back in port if needed
	if (nout > 0) {
		simLog(NULL,11,"ERROR: failure in nsv_integrate (shouldnt get here)\n");
		er = portputmols(sim,port,nout,species[0],species,positions);
		if (er != 0) simLog(NULL,11,"ERROR: failure in nsv_integrate (while putting particles back in the port)\n");
	}
}

vtkUnstructuredGrid* nsv_get_grid(NextSubvolumeMethod* nsv) {
	using namespace Kairos;

	return get_vtk_grid(nsv->get_grid(),nsv->get_diffusing_species());
}

