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

#define NSVC_CPP
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

	/*
	 * add any particles to lattice
	 */
	const int ns = nsv->get_diffusing_species().size();
	for (int var = 0; var < max; ++var) {

	}
	nsv->get_species(id)->particles.push_back(newr);

	//	const int ci = nsv->get_grid().get_cell_index(newr);
	//	nsv->get_species(id)->copy_numbers[ci]++;
	//	nsv->recalc_priority(ci);

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
		nsv->set_interface(interface,dt, false);

	} else if ((norm[1] == 1) || (norm[1] == -1)) {
		yrect interface(min,max,norm[1]);
		nsv->set_interface(interface,dt, false);

	} else if ((norm[2] == 1) || (norm[2] == -1)) {
		zrect interface(min,max,norm[2]);
		nsv->set_interface(interface,dt, false);

	}
}
void nsv_add_species(NextSubvolumeMethod* nsv,int id,double D) {
	using namespace Kairos;

	//TODO: assumes that species has not been added before
	Species *ns = new Species(id,D,&nsv->get_grid());
	nsv->add_diffusion(*ns);

}

class SmoldynSurface {
public:
	SmoldynSurface(surfaceptr s):s(s) {}
	bool lineXsurface(const Kairos::Vect3d& p1, const Kairos::Vect3d& p2) const {
		double crsspt[3];
		double p1t[3],p2t[3];
		p1t[0] = p1[0];
		p1t[1] = p1[1];
		p1t[2] = p1[2];
		p2t[0] = p2[0];
		p2t[1] = p2[1];
		p2t[2] = p2[2];
		for (int i = 0; i < PSMAX; ++i) {
			for (int j = 0; j < s->npanel[i]; ++j) {

				if (lineXpanel(p1t,p2t,s->panels[i][j],s->srfss->sim->dim,crsspt,NULL,NULL,NULL,NULL,NULL)==1) {
					return true;
				}
			}
		}
		return false;
	}
private:
	surfaceptr s;
};

class SmoldynCompartment {
public:
	SmoldynCompartment(compartptr c):c(c) {}
	bool is_in(const Kairos::Vect3d& test_point) const {
		double pt[3];
		pt[0] = test_point[0];
		pt[1] = test_point[1];
		pt[2] = test_point[2];
		return posincompart(c->cmptss->sim,pt,c)==1;
	}
private:
	compartptr c;
};

void nsv_add_reaction(NextSubvolumeMethod* nsv,rxnstruct *reaction) {
	const double rate = reaction->rate;
	const int nreactants = reaction->rxnss->order;
	const int *reactant_ids = reaction->rctident;
	const int nproducts = reaction->nprod;
	const int *product_ids = reaction->prdident;

	using namespace Kairos;

	ReactionSide lhs;
	for (int i = 0; i < nreactants; ++i) {
		lhs = lhs + *(nsv->get_species(reactant_ids[i]));
	}
	ReactionSide rhs;
	for (int i = 0; i < nproducts; ++i) {
		rhs = rhs + *(nsv->get_species(product_ids[i]));
	}

	if (reaction->srf) {
		SmoldynSurface surface(reaction->srf);
		nsv->add_reaction_on(rate, lhs >> rhs, surface);
	} else if (reaction->cmpt) {
		SmoldynCompartment compartment(reaction->cmpt);
		nsv->add_reaction_in(rate, lhs >> rhs, compartment);
	} else {
		nsv->add_reaction(rate, lhs >> rhs);
	}
}


//void nsv_add_reaction(NextSubvolumeMethod* nsv,double rate,
//						int nreactants,int *reactant_ids,
//						int nproducts,int *product_ids) {
//
//	//TODO: assumes that species has already been added via nsv_add_species
//	using namespace Kairos;
//
//
//	ReactionSide lhs;
//	for (int i = 0; i < nreactants; ++i) {
//		lhs = lhs + *(nsv->get_species(reactant_ids[i]));
//	}
//	ReactionSide rhs;
//	for (int i = 0; i < nproducts; ++i) {
//		rhs = rhs + *(nsv->get_species(product_ids[i]));
//	}
//	nsv->add_reaction(rate, lhs >> rhs);
//}

void nsv_integrate(NextSubvolumeMethod* nsv,double dt, portstruct *port) {
	using namespace Kairos;
	//std::cout << "running lattice dt"<<std::endl;
	nsv->integrate(dt);
	//std::cout << " lattice t = "<<nsv->get_time()<<" smoldyn time = "<<port->portss->sim->time<<std::endl;


	simptr sim = port->portss->sim;

	const int ns = nsv->get_diffusing_species().size();
	int n = 0;
	for (int i = 0; i < ns; ++i) {
		n += nsv->get_diffusing_species()[i]->particles.size();
	}
	n += sim->mols->nl[port->llport];

	int *species = new int[n];
	double **positions = new double*[n];


	/*
	 * look in port for new particles
	 */
	std::set<int> dirty_indicies;
	int nout = 0;
	n = sim->mols->nl[port->llport];
	for (int i = 0; i < n; ++i) {
		moleculeptr m = sim->mols->live[port->llport][i];

		//process other surface interactions
		int er = checksurfaces1mol(sim,m);
		if (er != 0) simLog(NULL,11,"ERROR: failure in nsv_integrate (while processing surface interactions)\n");

		//check if already gone back through port
		double *crsspt = new double[sim->dim];
		enum PanelFace face1,face2;
		bool in = true;
		for (int j = 0; j < port->srf->npanel[PSrect]; ++j) {
			if (lineXpanel(m->via,m->pos,port->srf->panels[PSrect][j],sim->dim,crsspt,&face1,&face2,NULL,NULL,NULL)) {
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

		// if outside lattice domain raise error
		Vect3d newr(0.5,0.5,0.5);
		for (int d = 0; d < sim->dim; ++d) {
			double low = nsv->get_grid().get_low()[d];
			double high = nsv->get_grid().get_high()[d];
			if (m->pos[d] < low) {
				simLog(NULL,11,"ERROR: particle unexpectedly outside lattice domain\n");
			} else if (m->pos[d] > high) {
				simLog(NULL,11,"ERROR: particle unexpectedly outside lattice domain\n");
			} else {
				newr[d] = m->pos[d];
			}
		}

		const int ci = nsv->get_grid().get_cell_index(newr);
		nsv->get_species(m->ident)->copy_numbers[ci]++;
		dirty_indicies.insert(ci);
	}
	for (std::set<int>::iterator i=dirty_indicies.begin();i!=dirty_indicies.end();i++) {
		nsv->recalc_priority(*i);
	}


	//delete particles in port
	int er = portgetmols(sim,port,-1,MSall,1);
	if (er != n) simLog(NULL,11,"ERROR: failure in nsv_integrate (while deleting particles from port)\n");


	/*
	 * put new particles in port
	 */
	for (int i = 0; i < ns; ++i) {
		Species *s = nsv->get_diffusing_species()[i];
		const int np = s->particles.size();
		for (int j = 0; j < np; ++j) {
			species[nout] = s->id;
			positions[nout] = s->particles[j].data();
			nout++;
		}
		s->particles.clear();
	}

	//put particles back in port if needed
	if (nout > 0) {
		er = portputmols(sim,port,nout,species[0],species,positions);
		if (er != 0) simLog(NULL,11,"ERROR: failure in nsv_integrate (while putting particles in the port)\n");
	}

	delete species;
	delete positions;
}

vtkUnstructuredGrid* nsv_get_grid(NextSubvolumeMethod* nsv) {
	using namespace Kairos;

	return get_vtk_grid(nsv->get_grid(),nsv->get_diffusing_species());
}

void nsv_molcountspace(NextSubvolumeMethod* nsv,int id, double *low, double *high, int dim, int nbins, int axis, int *ret_array) {
	using namespace Kairos;

	Vect3d vlow(0,0,0);
	Vect3d vhigh(1,1,1);
	Vect3d grid_size(1,1,1);

	for (int i = 0; i < dim; ++i) {
		vlow[i] = low[i];
		vhigh[i] = high[i];
		grid_size[i] = high[i] - low[i];
	}
	if (nbins > 1) {
		grid_size[axis] = (high[axis] - low[axis])/nbins;
	}

	StructuredGrid grid(vlow,vhigh,grid_size);

	std::vector<double> concentration;

	nsv->get_species(id)->get_concentration(grid,concentration);

	for (int i = 0; i < nbins; ++i) {
		ret_array[i] = concentration[i]*grid.get_cell_volume(i);
	}
}

void nsv_molcount(NextSubvolumeMethod* nsv, int *ret_array) {
	using namespace Kairos;

	std::vector<Species*> species = nsv->get_diffusing_species();
	for (unsigned int i = 0; i < species.size(); ++i) {
		const int n = std::accumulate(species[i]->copy_numbers.begin(),species[i]->copy_numbers.end(),0) +
						species[i]->particles.size();
		ret_array[species[i]->id] = n;
	}
}

void nsv_add_mol(NextSubvolumeMethod* nsv,int id, double* pos, int dim) {
	using namespace Kairos;
	Vect3d newr(0.5,0.5,0.5);
	// if outside lattice domain raise error
	for (int d = 0; d < dim; ++d) {
		double low = nsv->get_grid().get_low()[d];
		double high = nsv->get_grid().get_high()[d];
		if (pos[d] < low) {
			simLog(NULL,11,"ERROR: particle unexpectedly outside lattice domain\n");
		} else if (pos[d] > high) {
			simLog(NULL,11,"ERROR: particle unexpectedly outside lattice domain\n");
		} else {
			newr[d] = pos[d];
		}
	}
	//nsv->get_species(id)->particles.push_back(newr);

	const int ci = nsv->get_grid().get_cell_index(newr);
	nsv->get_species(id)->copy_numbers[ci]++;
	nsv->recalc_priority(ci);
}
