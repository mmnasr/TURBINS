/* Just a junk velocity field */
void Debugger_velocity_field(Velocity *vel, MAC_grid *grid) {
	
	int L_xs, L_ys, L_zs;
	int L_xe, L_ye, L_ze;
	int i, j, k;
	double ***data;
	int ierr; 
	
	/* Index of left-bottom-back corner on current processor for the DA data layout */
	L_xs = grid->L_xs;
	L_ys = grid->L_ys;
	L_zs = grid->L_zs;

	/* Index of right-top-front corner on current processor for the DA data layout */
	L_xe = grid->L_xe;
	L_ye = grid->L_ye;
	L_ze = grid->L_ze;

	ierr = DAVecGetArray(grid->DA_3D, vel->G_data, (void ***)&data); CHKERRQ(ierr);
	for (k=L_zs; k<L_ze; k++) {
		for (j=L_ys; j<L_ye; j++) {
			for (i=L_xs; i<L_xe; i++) {
				
				data[k][j][i] = (double)i*j*k;
					
			}
		}
	}

	ierr = DAVecRestoreArray(grid->DA_3D, vel->G_data, (void ***)&data); CHKERRQ(ierr);
	
}

	/* for debugging purposes */
void GVG_solve_Poisson_equation(Pressure *p, MAC_grid *grid, Parameters *params) {
	
	int L_xs, L_ys, L_zs;
	int L_xe, L_ye, L_ze;
	int i, j, k;
	double ***rhs;
	int ierr; 
	int index;
	int low, high; 
	int cell_lsys_index;
	double rhs_value;
	int NX, NY, NZ;
	PetscViewer viewer;
	
	NX = p->G_NX;
	NY = p->G_NY;
	NZ = p->G_NZ;
	/* Index of left-bottom-back corner on current processor for the DA data layout */
	L_xs = grid->L_xs;
	L_ys = grid->L_ys;
	L_zs = grid->L_zs;

	/* Index of right-top-front corner on current processor for the DA data layout */
	L_xe = grid->L_xe;
	L_ye = grid->L_ye;
	L_ze = grid->L_ze;

	PetscViewerASCIIOpen(PCW, "pResults.txt", &viewer);
	ierr = DAVecGetArray(grid->DA_3D, p->G_b, (void ***)&rhs); CHKERRQ(ierr);
	index = 0;
	
	ierr = VecGetOwnershipRange(p->G_b, &low, &high); CHKERRQ(ierr);
	
	MatView(p->A, 0);

	PetscSynchronizedPrintf(PCW, "Rank:%d Index_Low:%d Index_high:%d\n", params->rank, low, high);
	PetscSynchronizedFlush(PCW);
	
	for (k=L_zs; k<L_ze; k++) {
		for (j=L_ys; j<L_ye; j++) {
			for (i=L_xs; i<L_xe; i++) {
				
				cell_lsys_index = Pressure_get_lsys_index(p, i, j, k, 'p');
				if (!Grid_c_is_solid(grid, i, j, k)) {
				
					//rhs_value = (double) (k*NX*NY + j*NX + i);
					rhs[k][j][i] = (double) i*j*k + 1.0;
				}
				else {
					rhs[k][j][i] = 0.0;
				}
				//rhs[k][j][i] = rhs_value;

				//ierr = VecSetValue(p->G_b, cell_lsys_index, rhs_value, INSERT_VALUES);
			}
		}
	}

	//ierr = VecAssemblyBegin(p->G_b); CHKERRQ(ierr);
	//ierr = VecAssemblyEnd(p->G_b); CHKERRQ(ierr);

	ierr = DAVecRestoreArray(grid->DA_3D, p->G_b, (void ***)&rhs); CHKERRQ(ierr);
	
	//GVG_update_proc_ghost_nodes(&grid->DA_3D, &p->G_b, &p->L_b, 'I');
	VecCopy(p->G_b, p->G_data);
	
	int nlocal;
	double *array;
  	VecGetLocalSize(p->G_b, &nlocal);
   	VecGetArray(p->G_b,&array);
 	for (i=0; i<nlocal; i++) {
 		
		//printf("Rank:%d i:%d v:%f\n", params->rank, i, array[i]);
	}
	VecRestoreArray(p->G_b,&array);

	//Display_DA_3D_data(p->G_data, grid, params, "pressure", 'p');

	Pressure_solve(p);

/*
	double norm1A, norm2A, normIA;
	double norm1b, norm2b, normIb;
	double norm1x, norm2x, normIx;

	MatNorm(p->A, NORM_1, &norm1A);
	MatNorm(p->A, NORM_FROBENIUS, &norm2A);
	MatNorm(p->A, NORM_INFINITY, &normIA);

	VecNorm(p->G_b, NORM_1, &norm1b);
	VecNorm(p->G_b, NORM_2, &norm2b);
	VecNorm(p->G_b, NORM_INFINITY, &normIb);

	VecNorm(p->G_data, NORM_1, &norm1x);
	VecNorm(p->G_data, NORM_2, &norm2x);
	VecNorm(p->G_data, NORM_INFINITY, &normIx);
	
	PetscPrintf(PCW, "A  : norm1:%f norm2:%f normI:%f\n", norm1A, norm2A, normIA);
	PetscPrintf(PCW, "rhs: norm1:%f norm2:%f normI:%f\n", norm1b, norm2b, normIb);
	PetscPrintf(PCW, "sol: norm1:%f norm2:%f normI:%f\n", norm1x, norm2x, normIx);



	//Display_DA_3D_data(p->G_data, grid, params, "pressure", 'p');
	

	//MatView(p->A, 0);
	//VecView(p->G_b, 0);
	
	//Display_DA_3D_data(p->G_data, grid, params, "pressure", 'p');
	//getchar();
*/	

	PetscViewerASCIIPrintf(viewer, "Solution\n");
	VecView(p->G_data, viewer);
	PetscViewerDestroy(viewer);
}

	/* for debugging purposes */
void GVG_solve_vel_Poisson_equation(Velocity *vel, MAC_grid *grid, Parameters *params) {
	
	int L_xs, L_ys, L_zs;
	int L_xe, L_ye, L_ze;
	int i, j, k;
	double ***rhs;
	int ierr; 
	int index;
	int low, high; 
	int cell_lsys_index;
	double rhs_value;
	int NX, NY, NZ;
	short int (*velocity_is_solid)(MAC_grid *, int, int, int);
	PetscViewer viewer;
	
	switch (vel->component) {
		
		case 'u': velocity_is_solid = Grid_u_is_solid; 	PetscViewerASCIIOpen(PCW, "uResults.txt", &viewer); break;
		case 'v': velocity_is_solid = Grid_v_is_solid; 	PetscViewerASCIIOpen(PCW, "vResults.txt", &viewer); break;
		case 'w': velocity_is_solid = Grid_w_is_solid; 	PetscViewerASCIIOpen(PCW, "wResults.txt", &viewer); break;
	
	}
		
	NX = vel->G_NX;
	NY = vel->G_NY;
	NZ = vel->G_NZ;
	/* Index of left-bottom-back corner on current processor for the DA data layout */
	L_xs = grid->L_xs;
	L_ys = grid->L_ys;
	L_zs = grid->L_zs;

	/* Index of right-top-front corner on current processor for the DA data layout */
	L_xe = grid->L_xe;
	L_ye = grid->L_ye;
	L_ze = grid->L_ze;

	ierr = DAVecGetArray(grid->DA_3D, vel->G_b, (void ***)&rhs); CHKERRQ(ierr);
	index = 0;
	
	ierr = VecGetOwnershipRange(vel->G_b, &low, &high); CHKERRQ(ierr);
	
	MatView(vel->A, 0);

	PetscSynchronizedPrintf(PCW, "Rank:%d Index_Low:%d Index_high:%d\n", params->rank, low, high);
	PetscSynchronizedFlush(PCW);
	
	for (k=L_zs; k<L_ze; k++) {
		for (j=L_ys; j<L_ye; j++) {
			for (i=L_xs; i<L_xe; i++) {
				
				if (!velocity_is_solid(grid, i, j, k)) {
				
					//rhs_value = (double) (k*NX*NY + j*NX + i);
					rhs[k][j][i] = (double) i*j*k + 1.0;
				}
				else {
					rhs[k][j][i] = 0.0;
				}

			}
		}
	}

	ierr = DAVecRestoreArray(grid->DA_3D, vel->G_b, (void ***)&rhs); CHKERRQ(ierr);
	
	//VecCopy(vel->G_b, vel->G_data);
	
	int nlocal;
	double *array;
  	VecGetLocalSize(vel->G_b, &nlocal);
   	VecGetArray(vel->G_b,&array);
 	for (i=0; i<nlocal; i++) {
 		
		//printf("Rank:%d i:%d v:%f\n", params->rank, i, array[i]);
	}
	VecRestoreArray(vel->G_b,&array);


	Velocity_solve(vel);
/*

	Display_DA_3D_data(vel->G_data, grid, params, "velocity", vel->component);

	double norm1A, norm2A, normIA;
	double norm1b, norm2b, normIb;
	double norm1x, norm2x, normIx;

	MatNorm(vel->A, NORM_1, &norm1A);
	MatNorm(vel->A, NORM_FROBENIUS, &norm2A);
	MatNorm(vel->A, NORM_INFINITY, &normIA);

	VecNorm(vel->G_b, NORM_1, &norm1b);
	VecNorm(vel->G_b, NORM_2, &norm2b);
	VecNorm(vel->G_b, NORM_INFINITY, &normIb);

	VecNorm(vel->G_data, NORM_1, &norm1x);
	VecNorm(vel->G_data, NORM_2, &norm2x);
	VecNorm(vel->G_data, NORM_INFINITY, &normIx);
	
	PetscPrintf(PCW, "A  : norm1:%f norm2:%f normI:%f\n", norm1A, norm2A, normIA);
	PetscPrintf(PCW, "rhs: norm1:%f norm2:%f normI:%f\n", norm1b, norm2b, normIb);
	PetscPrintf(PCW, "sol: norm1:%f norm2:%f normI:%f\n", norm1x, norm2x, normIx);



	//Display_DA_3D_data(p->G_data, grid, params, "pressure", 'p');
	

	//MatView(p->A, 0);
	//VecView(p->G_b, 0);
	
	//Display_DA_3D_data(p->G_data, grid, params, "pressure", 'p');
	//getchar();
*/


	PetscViewerASCIIPrintf(viewer, "Solution\n");
	VecView(vel->G_data, viewer);
	PetscViewerDestroy(viewer);
	
}

