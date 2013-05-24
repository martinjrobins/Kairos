/* 
 * nsvc.h
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

#ifndef NSVC_H_
#define NSVC_H_

struct portstruct;

#ifdef __cplusplus
	#ifdef HAVE_VTK
		#include <vtkUnstructuredGrid.h>
		#include <string>
		#include "Kairos.h"
		using Kairos::NextSubvolumeMethod;
	#else
		typedef struct vtkUnstructuredGrid vtkUnstructuredGrid;
		typedef struct NextSubvolumeMethod NextSubvolumeMethod;
	#endif
#else
	typedef struct vtkUnstructuredGrid vtkUnstructuredGrid;
	typedef struct NextSubvolumeMethod NextSubvolumeMethod;
#endif


#ifdef __cplusplus
extern "C" {
#endif

extern void nsv_init();

extern NextSubvolumeMethod* nsv_new(double *min, double *max, double *dx, int dim);
extern void nsv_delete(NextSubvolumeMethod* nsv);
extern void nsv_print(NextSubvolumeMethod* nsv, char *buffer);
extern void nsv_add_interface(NextSubvolumeMethod* nsv,double dt, double *start,double *end,double *norm,int dim);
extern void nsv_add_species(NextSubvolumeMethod* nsv,int id,double D);
extern void nsv_add_reaction(NextSubvolumeMethod* nsv,double rate,
						int nreactants,int *reactant_ids,
						int nproducts,int *product_ids);
extern void nsv_integrate(NextSubvolumeMethod* nsv,double dt, portstruct *port);
extern vtkUnstructuredGrid* nsv_get_grid(NextSubvolumeMethod* nsv);

#ifdef __cplusplus
}
#endif

#endif /* NSVC_H_ */
