module global_data
    integer,parameter :: RKD = 8



    real(kind=RKD),parameter :: PI = 4.0*atan(1.0)

    real(kind=RKD),parameter :: SMV = tiny(real(1.0,8)) 

    real(kind=RKD),parameter :: UP = 1.0 

    real(kind=RKD) :: cfl 

    real(kind=RKD) :: dt

    real(kind=RKD) :: sim_time 

    real(kind=RKD) :: max_time 

! case1   character(len=20),parameter :: HSTFILENAME = "Sod2.hst" 

! case1   character(len=20),parameter :: RSTFILENAME = "Sod2.dat" 
    
   character(len=20),parameter :: HSTFILENAME = "Sjog2.hst" 

   character(len=20),parameter :: RSTFILENAME = "Sjog2.dat" 
    
    integer :: iter 





    real(kind=RKD) :: gam 

    integer :: ck 




    integer,parameter :: HSTFILE = 20

    integer,parameter :: RSTFILE = 21 

 

    type :: cell_center



        real(kind=RKD) :: x 

        real(kind=RKD) :: length 
        real(kind=RKD) :: w(3)  
        real(kind=RKD) :: w_x(3)
        real(kind=RKD) :: prim(3)

    end type cell_center
    
    type :: interface_average
        real(kind=RKD):: inter_w(3)
        real(kind=RKD):: inter_prim(3)
        real(kind=RKD):: inter_a(3)
    end type interface_average



    type :: cell_interface

        real(kind=RKD) :: flux_av(3) 
        real(kind=RKD) :: flux_a(3)

    end type cell_interface




    integer :: ixmin,ixmax 

    type(cell_center),allocatable,dimension(:) :: ctr 

    type(cell_interface),allocatable,dimension(:) :: vface 
 
    type(interface_average), allocatable, dimension(:) ::vf_L, vf_R

end module global_data





module tools

    use global_data

    implicit none

    contains


        function get_conserved(prim)

            real(kind=RKD),intent(in) :: prim(3)

            real(kind=RKD) :: get_conserved(3)



            get_conserved(1) = prim(1)

            get_conserved(2) = prim(1)*prim(2)

            get_conserved(3) = prim(3)/(gam-1.0)+0.5*prim(1)*prim(2)**2

        end function get_conserved



    

        function get_primary(w)

            real(kind=RKD),intent(in) :: w(3)

            real(kind=RKD) :: get_primary(3) 



            get_primary(1) = w(1)

            get_primary(2) = w(2)/w(1)

            get_primary(3) = (gam-1.0)*(w(3)-0.5*w(2)**2/w(1))

        end function get_primary




       

        function get_gamma(ck)

            integer,intent(in) :: ck

            real(kind=RKD) :: get_gamma



            get_gamma = float(ck+3)/float(ck+1)

        end function get_gamma



        function get_sos(prim)

            real(kind=RKD),intent(in) :: prim(3)

            real(kind=RKD) :: get_sos 



            get_sos = sqrt(gam*prim(3)/prim(1))

        end function get_sos



end module tools




module flux

    use global_data

    use tools

    implicit none



    integer,parameter :: MNUM = 6 

    integer,parameter :: MTUM = 4 



    contains


        subroutine calc_flux(cell_L,face,cell_R)

            type(interface_average),intent(in) :: cell_L,cell_R

            type(cell_interface),intent(inout) :: face

            real(kind=RKD) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<\xi^l>


!!!face%flux(1) = cell_L%prim(1)*( cell_L%prim(2)/2.0*erfc(-sqrt(cell_L%prim(3))*cell_L%prim(2)) &

!!!                                +0.5*exp(-cell_L%prim(3)*cell_L%prim(2)**2)/sqrt(PI*cell_L%prim(3)) ) &

!!!              +cell_R%prim(1)*( cell_R%prim(2)/2.0*erfc(sqrt(cell_R%prim(3))*cell_R%prim(2)) &

!!!                                -0.5*exp(-cell_R%prim(3)*cell_R%prim(2)**2)/sqrt(PI*cell_R%prim(3)) )

!!!

!!!face%flux(2) = cell_L%prim(1)*( (cell_L%prim(2)**2/2.0+1.0/4.0/cell_L%prim(3))*erfc(-sqrt(cell_L%prim(3))*cell_L%prim(2)) &

!!!                                +0.5*cell_L%prim(2)*exp(-cell_L%prim(3)*cell_L%prim(2)**2)/sqrt(PI*cell_L%prim(3)) ) &

!!!              +cell_R%prim(1)*( (cell_R%prim(2)**2/2.0+1.0/4.0/cell_R%prim(3))*erfc(sqrt(cell_R%prim(3))*cell_R%prim(2)) &

!!!                                -0.5*cell_R%prim(2)*exp(-cell_R%prim(3)*cell_R%prim(2)**2)/sqrt(PI*cell_R%prim(3)) )

!!!

!!!face%flux(3) = cell_L%prim(1)*( (cell_L%prim(2)**3/4.0+(ck+3)/8.0/cell_L%prim(3)*cell_L%prim(2)) &

!!!                                  *erfc(-sqrt(cell_L%prim(3))*cell_L%prim(2)) &

!!!                               +(cell_L%prim(2)**2/4.0+(ck+2)/8.0/cell_L%prim(3)) &

!!!                                  *exp(-cell_L%prim(3)*cell_L%prim(2)**2)/sqrt(PI*cell_L%prim(3)) ) &

!!!              +cell_R%prim(1)*( (cell_R%prim(2)**3/4.0+(ck+3)/8.0/cell_R%prim(3)*cell_R%prim(2)) &

!!!                                  *erfc(sqrt(cell_R%prim(3))*cell_R%prim(2)) &

!!!                                -(cell_R%prim(2)**2/4.0+(ck+2)/8.0/cell_R%prim(3)) &

!!!                                  *exp(-cell_R%prim(3)*cell_R%prim(2)**2)/sqrt(PI*cell_R%prim(3)) )
        




        call moment_u(cell_L%inter_prim,Mu,Mxi,Mu_L,Mu_R)

        face%flux_av(1) = cell_L%inter_prim(1)*Mu_L(1)

        face%flux_av(2) = cell_L%inter_prim(1)*Mu_L(2)

        face%flux_av(3) = 0.5*cell_L%inter_prim(1)*(Mu_L(3)+Mu_L(1)*Mxi(1))

        
        face%flux_a(1) =cell_L%inter_prim(1)*(cell_L%inter_a(1)*Mu_L(2)+cell_L%inter_a(2)*Mu_L(3)+0.5*cell_L%inter_a(3)*(Mu_L(4)+Mu_L(2)*Mxi(1)))
        
        face%flux_a(2) =cell_L%inter_prim(1)*(cell_L%inter_a(1)*Mu_L(3)+cell_L%inter_a(2)*Mu_L(4)+0.5*cell_L%inter_a(3)*(Mu_L(5)+Mu_L(3)*Mxi(1)))
        
       face%flux_a(3) =0.5*cell_L%inter_prim(1)*(cell_L%inter_a(1)*(Mu_L(4)+Mu_L(2)*Mxi(1))+cell_L%inter_a(2)*(Mu_L(5)+Mu_L(3)*Mxi(1))+0.5*cell_L%inter_a(3)*(Mu_L(6)+2*Mu_L(4)*Mxi(1)+Mu_L(2)*Mxi(2)))
        
        
        



        call moment_u(cell_R%inter_prim,Mu,Mxi,Mu_L,Mu_R)

        face%flux_av(1) = face%flux_av(1)+cell_R%inter_prim(1)*Mu_R(1)

        face%flux_av(2) = face%flux_av(2)+cell_R%inter_prim(1)*Mu_R(2)

        face%flux_av(3) = face%flux_av(3)+0.5*cell_R%inter_prim(1)*(Mu_R(3)+Mu_R(1)*Mxi(1))
    
        
        face%flux_a(1) = face%flux_a(1)+cell_R%inter_prim(1)*(cell_R%inter_a(1)*Mu_R(2)+cell_R%inter_a(2)*Mu_R(3)+0.5*cell_R%inter_a(3)*(Mu_R(4)+Mu_R(2)*Mxi(1)))
        
        face%flux_a(2) =face%flux_a(2)+cell_R%inter_prim(1)*(cell_R%inter_a(1)*Mu_R(3)+cell_R%inter_a(2)*Mu_R(4)+0.5*cell_R%inter_a(3)*(Mu_R(5)+Mu_R(3)*Mxi(1)))
        
        face%flux_a(3) =face%flux_a(3)+0.5*cell_R%inter_prim(1)*(cell_R%inter_a(1)*(Mu_R(4)+Mu_R(2)*Mxi(1))+cell_R%inter_a(2)*(Mu_R(5)+Mu_R(3)*Mxi(1))+0.5*cell_R%inter_a(3)*(Mu_R(6)+2*Mu_R(4)*Mxi(1)+Mu_R(2)*Mxi(2)))
         





        end subroutine calc_flux






        subroutine moment_u(prim,Mu,Mxi,Mu_L,Mu_R)

            real(kind=RKD),intent(in) :: prim(3)

            real(kind=RKD),intent(out) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM)

            real(kind=RKD),intent(out) :: Mxi(0:2)

            integer :: i



          

            Mu_L(0) = 0.5*erfc(-sqrt(0.5*prim(1)/prim(3))*prim(2))

            Mu_L(1) = prim(2)*Mu_L(0)+0.5*exp(-0.5*prim(1)/prim(3)*prim(2)**2)/sqrt(PI*0.5*prim(1)/prim(3))

            Mu_R(0) = 0.5*erfc(sqrt(0.5*prim(1)/prim(3))*prim(2))

            Mu_R(1) = prim(2)*Mu_R(0)-0.5*exp(-0.5*prim(1)/prim(3)*prim(2)**2)/sqrt(PI*0.5*prim(1)/prim(3))



            do i=2,MNUM

                Mu_L(i) = prim(2)*Mu_L(i-1)+(i-1)*Mu_L(i-2)*prim(3)/prim(1)

                Mu_R(i) = prim(2)*Mu_R(i-1)+(i-1)*Mu_R(i-2)*prim(3)/prim(1)

            end do



            Mu = Mu_L+Mu_R



            

            Mxi(0) = 1.0 !<\xi^0>

            Mxi(1) = ck*prim(3)/prim(1) !<\xi^2>

            Mxi(2) = (ck**2.0+2.0*ck)*(prim(3)/prim(1))**2.0 

        end subroutine moment_u

end module flux





module solver

    use global_data

    use tools

    use flux

    implicit none

    contains


        subroutine timestep()

            real(kind=RKD) :: tmax

            real(kind=RKD) :: sos 

            real(kind=RKD) :: prim(3) 

            integer :: i



            

            tmax = 0.0



            do i=ixmin+1,ixmax+1

               

                prim = get_primary(ctr(i)%w)



               

                sos = get_sos(prim)



                

                prim(2) = Abs(prim(2))+sos



                

                tmax = max(tmax,prim(2)/ctr(i)%length)

            end do



    

            dt = cfl/tmax

        end subroutine timestep



       

        subroutine evolution()

            integer :: i



            do i=ixmin,ixmax+1
                call calc_flux(vf_L(i),vface(i),vf_R(i))

            end do



        end subroutine evolution



      

        subroutine update()

            integer :: i,j



            do i=ixmin+1,ixmax+1

             

                ctr(i)%w = ctr(i)%w+dt*((vface(i-1)%flux_av-vface(i)%flux_av)-dt*(vface(i-1)%flux_a-vface(i)%flux_a)/2.0)/ctr(i)%length

                ctr(i)%prim = get_primary(ctr(i)%w)

            end do
            
           ctr(ixmin)=ctr(ixmin+1)
            
            ctr(ixmin-1)=ctr(ixmin)
            
            ctr(ixmax+2)=ctr(ixmax+1)
            ctr(ixmax+3)=ctr(ixmax+2) 
            
            
        
       forall(j=1:3)
           forall(i=ixmin:ixmax+2, (ctr(i+1)%w(j)-ctr(i)%w(j))*(ctr(i)%w(j)-ctr(i-1)%w(j))<=0)
                     ctr(i)%w_x(j)=0.0d0
           endforall
           forall(i=ixmin:ixmax+2, (ctr(i+1)%w(j)-ctr(i)%w(j))*(ctr(i)%w(j)-ctr(i-1)%w(j))>0)
              ctr(i)%w_x(j)=2.0*(ctr(i+1)%w(j)-ctr(i)%w(j))*(ctr(i)%w(j)-ctr(i-1)%w(j))/ctr(i)%length/(ctr(i+1)%w(j)-ctr(i-1)%w(j))
          endforall
           endforall

            do i=ixmin,ixmax+1
            vf_L(i)%inter_w=ctr(i)%w+0.5*ctr(i)%length*ctr(i)%w_x
            vf_L(i)%inter_prim=get_primary(vf_L(i)%inter_w)
!            if (vf_L(i)%inter_prim(3)*vf_L(i)%inter_prim(1)<0.0) then
 !           vf_L(i)%inter_w=ctr(i)%w
 !           vf_L(i)%inter_prim=ctr(i)%prim
  !          end if 
            
            vf_R(i)%inter_w=ctr(i+1)%w-0.5*ctr(i+1)%length*ctr(i+1)%w_x
            vf_R(i)%inter_prim=get_primary(vf_R(i)%inter_w)
   !          if (vf_R(i)%inter_prim(3)*vf_R(i)%inter_prim(1)<0.0) then
    !        vf_R(i)%inter_w=ctr(i+1)%w
     !       vf_R(i)%inter_prim=ctr(i+1)%prim
      !      end if
            
            
            end do
          
            forall(i=ixmin:ixmax+1)
            
                vf_L(i)%inter_a(3)=vf_L(i)%inter_w(1)*(2.0*ctr(i)%w_x(3)+vf_L(i)%inter_prim(2)**2.0*ctr(i)%w_x(1)-2.0*vf_L(i)%inter_prim(2)*ctr(i)%w_x(2))/(ck+1.0)*(vf_L(i)%inter_prim(3)**2.0)-ctr(i)%w_x(1)/vf_L(i)%inter_prim(3)
                vf_L(i)%inter_a(2)=(ctr(i)%w_x(2)-vf_L(i)%inter_prim(2)*ctr(i)%w_x(1))/vf_L(i)%inter_prim(3)-vf_L(i)%inter_a(3)*vf_L(i)%inter_prim(2)
                vf_L(i)%inter_a(1)=ctr(i)%w_x(1)/vf_L(i)%inter_w(1)-vf_L(i)%inter_a(2)*vf_L(i)%inter_prim(2)-vf_L(i)%inter_a(3)*(vf_L(i)%inter_prim(2)**2.0/2.0+(ck+1.0)*vf_L(i)%inter_prim(3)/2.0/vf_L(i)%inter_w(1))
            endforall
            forall(i=ixmin:ixmax+1)
            vf_R(i)%inter_a(3)=vf_R(i)%inter_w(1)*(2.0*ctr(i+1)%w_x(3)+vf_R(i)%inter_prim(2)**2.0*ctr(i+1)%w_x(1)-2.0*vf_R(i)%inter_prim(2)*ctr(i+1)%w_x(2))/(ck+1.0)*(vf_R(i)%inter_prim(3)**2.0)-ctr(i+1)%w_x(1)/vf_R(i)%inter_prim(3)
            vf_R(i)%inter_a(2)=(ctr(i+1)%w_x(2)-vf_R(i)%inter_prim(2)*ctr(i+1)%w_x(1))/vf_R(i)%inter_prim(3)-vf_R(i)%inter_a(3)*vf_R(i)%inter_prim(2)
            vf_R(i)%inter_a(1)=ctr(i+1)%w_x(1)/vf_R(i)%inter_w(1)-vf_R(i)%inter_a(2)*vf_R(i)%inter_prim(2)-vf_R(i)%inter_a(3)*(vf_R(i)%inter_prim(2)**2.0/2.0+(ck+1.0)*vf_R(i)%inter_prim(3)/2.0/vf_R(i)%inter_w(1))
        
            endforall

        end subroutine update

end module solver




module io

    use global_data

    use tools
    
    implicit none
    

    contains

      

        subroutine init()
            real(kind=RKD) :: xlength 

            real(kind=RKD) :: xscale



            

            cfl = 0.5

           max_time = 0.15d0 
!     max_time =0.1d0



            

            ck = 4 

            gam = get_gamma(ck)



            

            xlength = 1.0

            xscale = 0.01

            call init_geometry(xlength,xscale) 

            call init_flow_field() 

        end subroutine init




        subroutine init_geometry(xlength,xscale)

            real(kind=RKD),intent(inout) :: xlength

            real(kind=RKD),intent(in) :: xscale

            integer :: xnum 

            real(kind=RKD) :: dx

            integer :: i



            

            xnum = nint(xlength/xscale)

            xlength = xnum*xscale

            write(*,*) "xnum=",xnum

            write(*,*) "xlength=",xlength



            

            ixmin = 1

            ixmax = xnum



            

            allocate(ctr(ixmin-1:ixmax+3)) 

            allocate(vface(ixmin:ixmax+1)) 
            allocate(vf_L(ixmin:ixmax+1))
            allocate(vf_R(ixmin:ixmax+1))


            

            dx = xlength/(ixmax-ixmin+1)

            write(*,*) "dx=",dx




            forall(i=ixmin-1:ixmax+3)

                ctr(i)%x = (i-1.5)*dx

                ctr(i)%length = dx

            end forall

        end subroutine init_geometry




        subroutine init_flow_field()

            real(kind=RKD) :: xscale

            real(kind=RKD) :: prim_L(3), prim_R(3) 

            real(kind=RKD) :: w_L(3), w_R(3)

            integer :: i, j

        w_L(1) = 1.0d0
         w_L(2) = 0.0d0
          w_L(3) = 2.5d0


    

 w_R(1) = 0.125d0 

     w_R(2) = 0.0d0

     w_R(3) = 0.25d0
            
            
            
 !       w_L(1) = 1.0d0

  !        w_L(2) = -2.0d0

  !       w_L(3) = 3.0d0


    

   ! w_R(1) = 1.0d0 

    !  w_R(2) = 2.0d0

   !  w_R(3) = 3.0d0



            prim_L = get_primary(w_L)

            prim_R = get_primary(w_R)



            forall(i=ixmin-1:(ixmin+ixmax)/2)

                ctr(i)%prim = prim_L

                ctr(i)%w = w_L
                
            endforall
            
            forall(i=(ixmin+ixmax)/2+1:ixmax+3)

                ctr(i)%prim = prim_R

                ctr(i)%w = w_R

            endforall
            
        
         forall(j=1:3)
             forall(i=ixmin:ixmax+2, (ctr(i+1)%w(j)-ctr(i)%w(j))*(ctr(i)%w(j)-ctr(i-1)%w(j))<=0)
              ctr(i)%w_x(j)=0.0d0      
             endforall
            forall(i=ixmin:ixmax+2, (ctr(i+1)%w(j)-ctr(i)%w(j))*(ctr(i)%w(j)-ctr(i-1)%w(j))>0)             
             ctr(i)%w_x(j)=2.0*(ctr(i+1)%w(j)-ctr(i)%w(j))*(ctr(i)%w(j)-ctr(i-1)%w(j))/xscale/(ctr(i+1)%w(j)-ctr(i-1)%w(j))
             endforall
          endforall

            do i=ixmin,ixmax+1
            vf_L(i)%inter_w=ctr(i)%w+0.5*ctr(i)%length*ctr(i)%w_x
            vf_L(i)%inter_prim=get_primary(vf_L(i)%inter_w)
        !    if (vf_L(i)%inter_prim(3)*vf_L(i)%inter_prim(1)<0.0) then
         !   vf_L(i)%inter_w=ctr(i)%w
          !  vf_L(i)%inter_prim=ctr(i)%prim
           ! end if
            
            
            vf_R(i)%inter_w=ctr(i+1)%w-0.5*ctr(i+1)%length*ctr(i+1)%w_x
            vf_R(i)%inter_prim=get_primary(vf_R(i)%inter_w)
    !        if (vf_R(i)%inter_prim(3)*vf_R(i)%inter_prim(1)<0.0) then
     !       vf_R(i)%inter_w=ctr(i+1)%w
      !      vf_R(i)%inter_prim=ctr(i+1)%prim
       !     end if
            
            
            end do
          
            forall(i=ixmin:ixmax+1)
            
                vf_L(i)%inter_a(3)=vf_L(i)%inter_w(1)*(2.0*ctr(i)%w_x(3)+vf_L(i)%inter_prim(2)**2.0*ctr(i)%w_x(1)-2.0*vf_L(i)%inter_prim(2)*ctr(i)%w_x(2))/(ck+1.0)*(vf_L(i)%inter_prim(3)**2.0)-ctr(i)%w_x(1)/vf_L(i)%inter_prim(3)
                vf_L(i)%inter_a(2)=(ctr(i)%w_x(2)-vf_L(i)%inter_prim(2)*ctr(i)%w_x(1))/vf_L(i)%inter_prim(3)-vf_L(i)%inter_a(3)*vf_L(i)%inter_prim(2)
                vf_L(i)%inter_a(1)=ctr(i)%w_x(1)/vf_L(i)%inter_w(1)-vf_L(i)%inter_a(2)*vf_L(i)%inter_prim(2)-vf_L(i)%inter_a(3)*(vf_L(i)%inter_prim(2)**2.0/2.0+(ck+1.0)*vf_L(i)%inter_prim(3)/2.0/vf_L(i)%inter_w(1))
            endforall
            forall(i=ixmin:ixmax+1)
            vf_R(i)%inter_a(3)=vf_R(i)%inter_w(1)*(2.0*ctr(i+1)%w_x(3)+vf_R(i)%inter_prim(2)**2.0*ctr(i+1)%w_x(1)-2.0*vf_R(i)%inter_prim(2)*ctr(i+1)%w_x(2))/(ck+1.0)*(vf_R(i)%inter_prim(3)**2.0)-ctr(i+1)%w_x(1)/vf_R(i)%inter_prim(3)
            vf_R(i)%inter_a(2)=(ctr(i+1)%w_x(2)-vf_R(i)%inter_prim(2)*ctr(i+1)%w_x(1))/vf_R(i)%inter_prim(3)-vf_R(i)%inter_a(3)*vf_R(i)%inter_prim(2)
            vf_R(i)%inter_a(1)=ctr(i+1)%w_x(1)/vf_R(i)%inter_w(1)-vf_R(i)%inter_a(2)*vf_R(i)%inter_prim(2)-vf_R(i)%inter_a(3)*(vf_R(i)%inter_prim(2)**2.0/2.0+(ck+1.0)*vf_R(i)%inter_prim(3)/2.0/vf_R(i)%inter_w(1))
        
            endforall  

        
            

        end subroutine init_flow_field

        subroutine output()

            integer :: i

         

            open(unit=RSTFILE,file=RSTFILENAME,status="unknown",action="write")

            write(RSTFILE,*) "VARIABLES = Position, Density, Velocity, Pressure"

            do i=ixmin+1,ixmax+1

                write(RSTFILE,*) ctr(i)%x, ctr(i)%prim(1), ctr(i)%prim(2), ctr(i)%prim(3)

            enddo

            

            close(RSTFILE)

        end subroutine output

end module io





program main

    use global_data

    use solver

    use io

    implicit none



    

    call init()





    iter = 1 

    sim_time = 0.0 





    open(unit=HSTFILE,file=HSTFILENAME,status="replace",action="write") 

    write(HSTFILE,*) "VARIABLES = iter, sim_time, dt"




    do while(.true.)

        call timestep()

        if((sim_time+dt)>=max_time) dt=max_time-sim_time

        call evolution() 

        call update() 





        if (mod(iter,10)==0) then

            write(*,"(A18,I15,2E15.7)") "iter,sim_time,dt:",iter,sim_time,dt

            write(HSTFILE,"(I15,2E15.7)") iter,sim_time,dt

        end if



        

        if (sim_time>=max_time) exit



        iter = iter+1

        sim_time = sim_time+dt

    end do



    write(*,"(A18,I15,2E15.7)") "iter,sim_time,dt:",iter,sim_time,dt

    write(HSTFILE,"(I15,2E15.7)") iter,sim_time,dt



    

    close(HSTFILE)



    

    call output()

end program main