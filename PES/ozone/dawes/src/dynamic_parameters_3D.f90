MODULE dynamic_parameters
implicit none
save
public
real*8, allocatable :: b2(:,:),b2_lower(:,:),b2_minimal(:,:),b2_seed(:,:),d(:),d_seed(:),Jac(:),Jac2(:),coords(:,:),coords_seed(:,:)
real*8, allocatable :: cart(:),bdist(:),stored_weights(:),stored_grad(:,:)
real*8, allocatable :: pot(:),pot_seed(:),grad(:,:),grad_seed(:,:),dcart(:)
integer :: basis_1,basis_2,basis_3,order_1,order_2,order_3,count3,zz,zz_low,zz4,support,count7,myid,natom,ab_flag,ab_flag2,lab,permfac
integer :: maxpoints,nbdist,order_1_min,order_2_min,order_3_min,count_seed,low_grid,subzero
real*8 :: epss,acc,rmax(3),rmin(3),poten,hartokcl,ugrad,Max_E,Max_E_seed,E_range,dist_tol(3),E_limit,alpha(3),Glob_min,current_geom(3),current_geom_grad(3),mass(3)
character(len=3) :: symb(3)
integer,allocatable :: stored_weights_ind(:)
END MODULE dynamic_parameters
