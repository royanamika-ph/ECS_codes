!!!For python wrapper run!!! $ python3 -m numpy.f2py -c -m aggregation aggregation.f90

subroutine pixelize(datas, n, matrix)
    real :: datas(:,:) ! data contains x, y and radius.
    integer :: n, m
    real :: L, x_b, y_b, x_p, y_p, radius, d
    integer(kind=1) :: matrix(n,n)
    !f2py real,intent(in) :: datas
    !f2py integer(kind=1), intent(out) :: matrix
    !f2py integer, intent(in) :: n
    
    L = 1.0
    l = L/n  ! pixel size = box length/no of pixels in a row
    m = size(datas(:,1)) ! number of elements
    iloop: do i=1,n
        jloop: do j =1,n
            x_b = l*(i-1) + l/2
            y_b = l*(j-1) + l/2
            matrix(j,i) = 1
            idloop: do id = 1,m
                x_p = datas(id,1)
                y_p = datas(id,2)
                radius = datas(id,3)
                d = sqrt((x_p-x_b)**2 + (y_p-y_b)**2)
                if (d < radius) then
                    matrix(j,i) = 0
                    exit idloop  !! if the pixel is within any one element, no need to check for other elements.
                endif
            enddo idloop
        enddo jloop
    enddo iloop
    
end subroutine

subroutine hexa_chi(matrix,n,chi_hex)
    implicit none
    !integer(kind=1)  :: matrix(:,:)
    integer,dimension(2*n,2*n) :: matrix
    integer  ::  chi_hex ,n 
  
    real,dimension(n,n) :: rhex
    integer,dimension(n,n):: label,sizeb,sizew,trialno
    integer,dimension(n,n)::check1, check2, check3, check4, check5, check6, bcheck1, bcheck2, pointvalue
    integer,dimension(n,n)::   wcheck1, wcheck2, wcheck3, wcheck4, wcheck5, wcheck6, xpoint, ypoint
    integer, dimension (n*n) :: iiw,jjw,ii,jj,nbn,nwn,clustersize,num,numf,trno,tr
    real,dimension(2) :: euler, eulerno,normeuler,u,p
    real :: i1,j1,j2,inc,nbf,nwf
    integer :: blacklabel ,nn ,nnw,sno,wsno ,totn ,totnw,sp,spw
    integer :: whitelabel,m1,m2,b1,b2,nb,nw,wl,euln,itrial,nbt,nwt
    integer :: tot,isz ,t1,trial
    integer :: i, j, k
          
    !f2py intent(in) :: matrix,n
    !f2py intent(out) :: chi_hex

    do k=1,1
        
        nbt=0
        nwt=0
        tot=1
          
        do isz=1,n*n
          num(isz)=0
        enddo
          
        trial=1
        t1=1
        do itrial =1,trial

          j2=0.0
          j1=0.0
          blacklabel=0
          whitelabel=0

          do i=1,n
            do j=1,n
              if(MOD(i,2).eq.0)then
                rhex(i,j)= matrix((2*(i-1)+1),(2*j))
              else
                rhex(i,j)=matrix((2*(i-1)+1),(2*(j-1)+1))
              endif
            enddo
          enddo

          do i=1,n
            do j=1,n
              if (rhex(i,j).eq.0)then
                pointvalue(i,j)=1
              else
                pointvalue(i,j)=0
              endif
            enddo
          enddo

          do i=1,n
            do j=1,n
              label(i,j)=0
            enddo
          enddo

          do i=2,n-1
            do j=2,n-1
              if ((pointvalue(i,j).eq.1 ).and. (label(i,j).eq.0))then

                blacklabel=blacklabel+1
                label(i,j)=blacklabel
                nn=0
                sno=0
                ii(sno)=i
                jj(sno)=j
            
                sp=0
11              totn=nn

                do sno=sp,totn
        
                  check1(ii(sno),jj(sno))=pointvalue(ii(sno),jj(sno)-1)
                  check2(ii(sno),jj(sno))=pointvalue(ii(sno),jj(sno)+1)
                  check3(ii(sno),jj(sno))=pointvalue(ii(sno)-1,jj(sno))
                  check4(ii(sno),jj(sno))=pointvalue(ii(sno)-1,jj(sno)+1)
                  check5(ii(sno),jj(sno))=pointvalue(ii(sno)+1,jj(sno)-1)
                  check6(ii(sno),jj(sno))=pointvalue(ii(sno)+1,jj(sno))
                  if (check1(ii(sno),jj(sno)).eq.1 )then
                    if(label(ii(sno),jj(sno)-1).eq.0)then
                      label(ii(sno),jj(sno)-1)=label(ii(sno),jj(sno))
                      if(ii(sno).gt.1 .and. ii(sno).lt.n )then
                        if((jj(sno)-1).gt.1 .and. (jj(sno)-1).lt.n) then
                          ii(nn+1)=ii(sno)
                          jj(nn+1)=jj(sno)-1
                          nn=nn+1
                        endif
                      endif
                    endif
                  endif
                  if(check2(ii(sno),jj(sno)).eq.1)then
                    if(label(ii(sno),jj(sno)+1).eq.0)then
                      label(ii(sno),jj(sno)+1)=label(ii(sno),jj(sno))
                      if(ii(sno).gt.1 .and. ii(sno).lt.n )then
                        if((jj(sno)+1).gt.1 .and. (jj(sno)+1).lt.n) then
                          ii(nn+1)=ii(sno)
                          jj(nn+1)=jj(sno)+1
                          nn=nn+1
                        endif
                      endif
                    endif
                  endif

                  if (check3(ii(sno),jj(sno)).eq.1 )then
                    if(label(ii(sno)-1,jj(sno)).eq.0)then
                      label(ii(sno)-1,jj(sno))=label(ii(sno),jj(sno))
                      if((ii(sno)-1).gt.1 .and. (ii(sno)-1).lt.n )then
                        if(jj(sno).gt.1 .and. jj(sno).lt.n) then
                          ii(nn+1)=ii(sno)-1
                          jj(nn+1)= jj(sno)
                          nn=nn+1
                        endif
                      endif
                    endif
                  endif

                  if (check4(ii(sno),jj(sno)).eq.1 )then
                    if(label(ii(sno)-1,jj(sno)+1).eq.0)then
                      label(ii(sno)-1,jj(sno)+1)=label(ii(sno),jj(sno))
                      if((ii(sno)-1).gt.1 .and. (ii(sno)-1).lt.n )then
                        if((jj(sno)+1).gt.1 .and. (jj(sno)+1).lt.n) then
                          ii(nn+1)=ii(sno)-1
                          jj(nn+1)=jj(sno)+1
                          nn=nn+1
                        endif
                      endif
                    endif
                  endif
                  
                  if (check5(ii(sno),jj(sno)).eq.1 )then
                    if(label(ii(sno)+1,jj(sno)-1).eq.0)then
                      label(ii(sno)+1,jj(sno)-1)=label(ii(sno),jj(sno))
                      if((ii(sno)+1).gt.1 .and. (ii(sno)+1).lt.n )then
                        if((jj(sno)-1).gt.1 .and.(jj(sno)-1).lt.n) then
                          ii(nn+1)=ii(sno)+1
                          jj(nn+1)=jj(sno)-1
                          nn=nn+1
                        endif
                      endif
                    endif
                  endif
                  
                  if (check6(ii(sno),jj(sno)).eq.1 )then
                    if(label(ii(sno)+1,jj(sno)).eq.0)then
                      label(ii(sno)+1,jj(sno))=label(ii(sno),jj(sno))
                      if((ii(sno)+1).gt.1 .and. (ii(sno)+1).lt.n )then
                        if(jj(sno).gt.1 .and. jj(sno).lt.n) then
                          ii(nn+1)=ii(sno)+1
                          jj(nn+1)=jj(sno)
                          nn=nn+1
                        endif
                      endif
                    endif
                  endif
 
 
                  if (sno.eq.0)then
                    nbn(sno)=nn
                  else
                    nbn(sno)=nn-nbn(sno-1)
                  endif
                  
                  if (nbn(sno).gt.0)then
                    sp=sp+1
                    goto 11
                  endif
                enddo
              endif
            enddo
          enddo



     
             do i=2,n-1
           do j=2,n-1
          ! write(*,*)i,j

           if ((pointvalue(i,j).eq.0 ).and. (label(i,j).eq.0))then

           whitelabel=whitelabel-1

           label(i,j)=whitelabel
           !write(*,*)'hi',point(i,j)
           !write(*,*)'label',label(i,j),point(i,j),i,j
           nnw=0
           wsno=0
           iiw(wsno)=i
           jjw(wsno)=j
           !write(*,*)ii(0),jj(0)
           spw=0


 33      totnw=nnw



           do wsno=spw,totnw
        ! if(iiw(wsno).gt.1 .and. iiw(wsno).lt.n )then
         ! if(jjw(wsno).gt.1 .and. jjw(wsno).lt.n) then



        wcheck1(iiw(wsno),jjw(wsno))=pointvalue(iiw(wsno),jjw(wsno)-1)
        wcheck2(iiw(wsno),jjw(wsno))=pointvalue(iiw(wsno),jjw(wsno)+1)
        wcheck3(iiw(wsno),jjw(wsno))=pointvalue(iiw(wsno)-1,jjw(wsno))
        wcheck4(iiw(wsno),jjw(wsno))=pointvalue(iiw(wsno)-1,jjw(wsno)+1)
        wcheck5(iiw(wsno),jjw(wsno))=pointvalue(iiw(wsno)+1,jjw(wsno)-1)
        wcheck6(iiw(wsno),jjw(wsno))=pointvalue(iiw(wsno)+1,jjw(wsno))
        if (wcheck1(iiw(wsno),jjw(wsno)).eq.0)then
        if(label(iiw(wsno),jjw(wsno)-1).eq.0)then
        label(iiw(wsno),jjw(wsno)-1)=label(iiw(wsno),jjw(wsno))

        if(iiw(wsno).gt.1 .and. iiw(wsno).lt.n )then
          if((jjw(wsno)-1).gt.1 .and. (jjw(wsno)-1).lt.n) then


        iiw(nnw+1)=iiw(wsno)
        jjw(nnw+1)=jjw(wsno)-1
         nnw=nnw+1


          ! goto 11
        endif
        endif
        endif
        endif
       if(wcheck2(iiw(wsno),jjw(wsno)).eq.0)then
       if(label(iiw(wsno),jjw(wsno)+1).eq.0)then
       label(iiw(wsno),jjw(wsno)+1)=label(iiw(wsno),jjw(wsno))

       if(iiw(wsno).gt.1 .and. iiw(wsno).lt.n )then
          if((jjw(wsno)+1).gt.1 .and. (jjw(wsno)+1).lt.n) then

        iiw(nnw+1)=iiw(wsno)
        jjw(nnw+1)=jjw(wsno)+1
        nnw=nnw+1

           !goto 11
        endif
        endif
        endif
        endif
        if (wcheck3(iiw(wsno),jjw(wsno)).eq.0 )then
        if(label(iiw(wsno)-1,jjw(wsno)).eq.0)then
        label(iiw(wsno)-1,jjw(wsno))=label(iiw(wsno),jjw(wsno))
        if((iiw(wsno)-1).gt.1 .and. (iiw(wsno)-1).lt.n )then
          if(jjw(wsno).gt.1 .and. jjw(wsno).lt.n) then

           iiw(nnw+1)=iiw(wsno)-1
           jjw(nnw+1)= jjw(wsno)
           nnw=nnw+1
           !write(*,*)ii(sno),jj(sno)
           !goto 11
        endif
        endif
        endif
        endif
        if (wcheck4(iiw(wsno),jjw(wsno)).eq.0 )then
        if(label(iiw(wsno)-1,jjw(wsno)+1).eq.0)then
           label(iiw(wsno)-1,jjw(wsno)+1)=label(iiw(wsno),jjw(wsno))
           if((iiw(wsno)-1).gt.1 .and. (iiw(wsno)-1).lt.n )then
          if((jjw(wsno)+1).gt.1 .and. (jjw(wsno)+1).lt.n) then

           iiw(nnw+1)=iiw(wsno)-1
           jjw(nnw+1)=jjw(wsno)+1
           nnw=nnw+1
           !goto 11
           endif
           endif
           endif
           endif
           if (wcheck5(iiw(wsno),jjw(wsno)).eq.0 )then
           if(label(iiw(wsno)+1,jjw(wsno)-1).eq.0)then
           label(iiw(wsno)+1,jjw(wsno)-1)=label(iiw(wsno),jjw(wsno))
           if((iiw(wsno)+1).gt.1 .and. (iiw(wsno)+1).lt.n )then
          if((jjw(wsno)-1).gt.1 .and. (jjw(wsno)-1).lt.n) then

           iiw(nnw+1)=iiw(wsno)+1
           jjw(nnw+1)=jjw(wsno)-1
           nnw=nnw+1
           !goto 11
           endif
           endif
           endif
           endif
           if (wcheck6(iiw(wsno),jjw(wsno)).eq.0 )then
           if(label(iiw(wsno)+1,jjw(wsno)).eq.0)then
           label(iiw(wsno)+1,jjw(wsno))=label(iiw(wsno),jjw(wsno))
           if((iiw(wsno)+1).gt.1 .and. (iiw(wsno)+1).lt.n )then
          if(jjw(wsno).gt.1 .and. jjw(wsno).lt.n) then

           iiw(nnw+1)=iiw(wsno)+1
           jjw(nnw+1)=jjw(wsno)
           nnw=nnw+1

           !goto 11
           endif
           endif
           endif
           endif



          if (wsno.eq.0)then
           nwn(wsno)=nnw
           !write(*,*)'jjj' ,nn,sno
           else
           nwn(wsno)=nnw-nwn(wsno-1)
          ! write(*,*),'gg',nwn(wsno)

           endif
          if (nwn(wsno).gt.0)then
         ! write(*,*)nbn(sno),sno,point(i,j)
          spw=spw+1
          !write(*,*)'kkk'
          goto 33
         ! write(*,*)'hi'

          endif


          !endif
          !endif


          enddo
          endif
          enddo
          enddo


! boundary points
             do i=1,1
             do j=1,n
             b1=0
             b2=0
              if(label(i,j).eq.0)then


              if (j.lt.n)then
              bcheck1(i,j)=pointvalue(i,j+1)
              endif
              if (j.gt.1) then
              bcheck2(i,j)=pointvalue(i,j-1)
              endif
              if (bcheck1(i,j).eq.pointvalue(i,j))then
              if(label(i,j+1).ne.0)then
              label(i,j)=label(i,j+1)
              b1=1
              endif
              endif
             if ( bcheck2(i,j).eq.pointvalue(i,j))then
             if(label(i,j-1).ne.0)then
             label(i,j)=label(i,j-1)
             b2=1
             endif
             endif
             if((b1.eq.1).and.(b2.eq.1)) then
             if(label(i,j-1).ne.label(i,j+1))then
             m1=min(label(i,j+1),label(i,j+1))
             label(i,j)=m1
             label(i,j-1)=m1
             label(i,j+1)=m1
             endif
             endif

             if (b1.eq.0 .and. b2.eq.0)then
             if (pointvalue(i,j).eq.1)then
             blacklabel=blacklabel+1
             label(i,j)=blacklabel

             else
             whitelabel=whitelabel-1
             label(i,j)=whitelabel
             endif
             endif
             endif
             enddo
             enddo


            do i=n,n
             do j=1,n
             b1=0
             b2=0
              if(label(i,j).eq.0)then


              if (j.lt.n)then
              bcheck1(i,j)=pointvalue(i,j+1)
              endif
              if (j.gt.1) then
              bcheck2(i,j)=pointvalue(i,j-1)
              endif
              if (bcheck1(i,j).eq.pointvalue(i,j))then
              if(label(i,j+1).ne.0)then
              label(i,j)=label(i,j+1)
              b1=1
              endif
              endif
             if ( bcheck2(i,j).eq.pointvalue(i,j))then
             if(label(i,j-1).ne.0)then
             label(i,j)=label(i,j-1)
             b2=1
             endif
             endif
             if((b1.eq.1).and.(b2.eq.1))then
             if(label(i,j-1).ne.label(i,j+1))then
             m1=min(label(i,j+1),label(i,j+1))
             label(i,j)=m1
             label(i,j-1)=m1
             label(i,j+1)=m1
             endif
             endif

             if (b1.eq.0 .and. b2.eq.0)then
             if (pointvalue(i,j).eq.1)then
             blacklabel=blacklabel+1
             label(i,j)=blacklabel

             else
             whitelabel=whitelabel-1
             label(i,j)=whitelabel
             endif
             endif
             endif
             enddo
             enddo



               do i=1,n
              j=1
             b1=0
             b2=0
              if(label(i,j).eq.0)then


              if (i.lt.n)then
              bcheck1(i,j)=pointvalue(i+1,j)
              endif
              if (i.gt.1) then
              bcheck2(i,j)=pointvalue(i-1,j)
              endif
              if (bcheck1(i,j).eq.pointvalue(i,j))then
              if(label(i+1,j).ne.0)then
              label(i,j)=label(i+1,j)
              b1=1
              endif
              endif
             if ( bcheck2(i,j).eq.pointvalue(i,j))then
             if(label(i-1,j).ne.0)then
             label(i,j)=label(i-1,j)
             b2=1
             endif
             endif
             if((b1.eq.1).and.(b2.eq.1))then
             if(label(i-1,j).ne.label(i+1,j))then
             m2=min(label(i-1,j),label(i+1,j))
             label(i,j)=m2
             label(i-1,j)=m2
             label(i+1,j)=m2
             endif
             endif

             if (b1.eq.0 .and. b2.eq.0)then
             if (pointvalue(i,j).eq.1)then
             blacklabel=blacklabel+1
             label(i,j)=blacklabel

             else
             whitelabel=whitelabel-1
             label(i,j)=whitelabel
             endif
             endif
             endif
             enddo



               do i=1,n
              j=n
             b1=0
             b2=0
              if(label(i,j).eq.0)then


              if (i.lt.n)then
              bcheck1(i,j)=pointvalue(i+1,j)
              endif
              if (i.gt.1) then
              bcheck2(i,j)=pointvalue(i-1,j)
              endif
              if (bcheck1(i,j).eq.pointvalue(i,j))then
              if(label(i+1,j).ne.0)then
              label(i,j)=label(i+1,j)
              b1=1
              endif
              endif
             if ( bcheck2(i,j).eq.pointvalue(i,j))then
             if(label(i-1,j).ne.0)then
             label(i,j)=label(i-1,j)
             b2=1
             endif
             endif
             if((b1.eq.1).and.(b2.eq.1))then
             if(label(i-1,j).ne.label(i+1,j))then
             m2=min(label(i-1,j),label(i+1,j))
             label(i,j)=m2
             label(i-1,j)=m2
             label(i+1,j)=m2
             endif
             endif

             if (b1.eq.0 .and. b2.eq.0)then
             if (pointvalue(i,j).eq.1)then
             blacklabel=blacklabel+1
             label(i,j)=blacklabel

             else
             whitelabel=whitelabel-1
             label(i,j)=whitelabel
             endif
             endif
             endif
             enddo
          nb=0
          wl=0
          do i=1,n
            do j=1,n
              if(label(i,j).gt.nb)then
                nb=nb+1
              endif
              if (label(i,j).lt.wl)then
                wl=wl-1
              endif
            enddo
          enddo

          nw = abs(wl)
          nbt=nbt+nb
          nwt=nwt+nw
        enddo  !trial ends
        nbf=nbt/trial
        nwf=nwt/trial
        euln=nbf-nwf
        chi_hex= euln
        write(*,*) nbf,nwf,euln,chi_hex
        euler(t1)=euln
        t1=t1+1
        eulerno(k)=euler(t1-1)
        u(k)=(eulerno(k))/(n*n)
        p(k+1)=p(k)+inc
    enddo
   
end subroutine hexa_chi

subroutine sq_chi(matrix,chi)
    !implicit doubleprecision(a-h,o-z)
    implicit none
    
    integer, parameter :: trial=1
    real, parameter :: prob=0.5
    integer(kind=1)  :: matrix(:,:)
    integer  ::  chi, lp
    integer, dimension(2*size(matrix,1), 2*size(matrix,1)) :: label
    real, dimension(size(matrix,1),size(matrix,1)) :: r2
    integer, dimension(size(matrix,1),size(matrix,1)) ::  latvalue, check1, check2, check3, check4, c1, c2, c3, c4
    integer, dimension(size(matrix,1)*size(matrix,1)) :: sizeb, sizew, clustersize, num, num1, tr ,trno, trialno 
    integer :: ab1,ab2,ab3,ab4,l1,l2,l3 ,i1,j1
    integer :: blacklabel ,voidlabel,m1,m2
    integer :: nb,nw,wl ,itrial ,nbt,nwt,nbf,nwf,eun,isz,tot,t1
    integer :: i, j, k, l, m, n
    real eun1
    
    !f2py intent(in) :: matrix
    !f2py intent(out) :: chi

    nbt=0
    nwt=0
    tot=1
    t1=1
    lp = size(matrix,1)

    do isz=1,lp*lp
        num(isz)=0
    enddo
    
    ! r == matrix
    do  itrial=1,trial
    
        do i=1,lp
            do j=1,lp
                blacklabel=0
                voidlabel=0
                check1(i,j)=-1
                check2(i,j)=-1
                check3(i,j)=-1
                check4(i,j)=-1
                c1(i,j)=0
                c2(i,j)=0
                c3(i,j)=0
                c4(i,j)=0

                if(matrix(i,j).eq.0)then
                    latvalue(i,j)=1
                    label(i,j)=blacklabel+1
                else
                    latvalue(i,j)=0
                    label(i,j)=voidlabel-1
                endif

                if(j.gt.1)then
                    check1(i,j)=latvalue(i,j-1)
                    if(check1(i,j).eq.latvalue(i,j))then
                        label(i,j)=label(i,j-1)
                        c1(i,j)=label(i,j-1)
                    endif
                endif

                if(i.gt.1)then
                    check2(i,j)=latvalue(i-1,j)
                    if(check2(i,j).eq.latvalue(i,j))then
                        label(i,j)=label(i-1,j)
                        c2(i,j)=label(i-1,j)
                    endif
                endif

                if((i.gt.1).and.(j.lt.lp))then
                    check3(i,j)=latvalue(i-1,j+1)
                    if((c2(i,j).eq.0).and.(c1(i,j).eq.0))then
                        if(check3(i,j).eq.latvalue(i,j)) then
                            r2(i,j)=rand( )
                            if(r2(i,j).ge.prob)then
                                label(i,j)= label(i-1,j+1)
                                c3(i,j)=label(i-1,j+1)
                            endif
                        endif
                    endif
                endif

                if((j.gt.1).and.(i.gt.1))then
                    check4(i,j)=latvalue(i-1,j-1)
                    if((c1(i,j).eq.0).and.(c2(i,j).eq.0))then
                        if(check4(i,j).eq.latvalue(i,j))then
                            if(c3(i,j-1).eq.0)then
                                label(i,j)=label(i-1,j-1)
                                c4(i,j)=label(i-1,j-1)
                            else
                                c4(i,j)=0
                            endif
                        endif
                    endif
                endif

                ab1=abs(c1(i,j))
                ab2=abs(c2(i,j))
                ab3=abs(c3(i,j))
                ab4=abs(c4(i,j))

                if((ab1.gt.0).and.(ab2.gt.0).and.(ab1.ne.ab2))then
                
                    l1=min(ab1,ab2)
                    i1=i
                    j1=j
                    if(latvalue(i,j).eq.0)then
                        label(i,j)=-l1
                    else
                        label(i,j)=l1
                    endif
                    
                    if (l1.eq.ab1)then
                    
                        do l=1,i1-1
                            do m=1,lp
                                if (label(l,m).eq.c2(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l1
                                    else
                                        label(l,m)=l1
                                    endif
                                endif
                            enddo
                        enddo

                        do l=i1,i1
                            do m=1,j1-1
                                if (label(l,m).eq.c2(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l1
                                    else
                                        label(l,m)=l1
                                    endif
                                endif
                            enddo
                        enddo
                        
                        if(latvalue(i,j).eq.0)then
                            do l=1,i1-1
                                do m=1,lp
                                    if ((label(l,m).ne.c1(i,j)).and.(label(l,m).ne.c2(i,j)))then
                                        if(label(l,m).lt.c2(i,j))then!modification
                                            label(l,m)=label(l,m)+1
                                        endif
                                    endif
                                enddo
                            enddo
                            do l=i1,i1
                                do m=1,j1-1
                                    if ((label(l,m).ne.c1(i,j)).and.(label(l,m).ne.c2(i,j)))then
                                        if(label(l,m).lt.c2(i,j))then !modification
                                            label(l,m)=label(l,m)+1
                                        endif
                                    endif
                                enddo
                            enddo
                        else
                            do l=1,i1-1
                                do m=1,lp
                                    if ((label(l,m).ne.c1(i,j)).and.(label(l,m).ne.c2(i,j)))then
                                        if(label(l,m).gt.c2(i,j))then
                                            label(l,m)=label(l,m)-1
                                        endif
                                    endif
                                enddo
                            enddo
                            do l=i1,i1
                                do m=1,j1-1
                                    if ((label(l,m).ne.c1(i,j)).and.(label(l,m).ne.c2(i,j)))then
                                        if(label(l,m).gt.c2(i,j))then
                                            label(l,m)=label(l,m)-1
                                        endif
                                    endif
                                enddo
                            enddo
                        endif

                    else
                    
                        do l=1,i1-1
                            do m=1,lp
                                if(label(l,m).eq.c1(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l1
                                    else
                                        label(l,m)=l1
                                    endif
                                endif
                            enddo
                        enddo

                        do l=i1,i1
                            do m=1,lp
                                if(label(l,m).eq.c1(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l1
                                    else
                                        label(l,m)=l1
                                    endif
                                endif
                            enddo
                        enddo

                        if(latvalue(i,j).eq.0)then
                            do l=1,i1-1
                                do m=1,lp
                                    if ((label(l,m).ne.c1(i,j)).and.(label(l,m).ne.c2(i,j)))then
                                        if(label(l,m).lt.c1(i,j))then!modification
                                            label(l,m)=label(l,m)+1
                                        endif
                                    endif
                                enddo
                            enddo
                            do l=i1,i1
                                do m=1,j1-1
                                    if ((label(l,m).ne.c1(i,j)) .and.(label(l,m).ne.c2(i,j)))then
                                        if(label(l,m).lt.c1(i,j))then!modification
                                            label(l,m)=label(l,m)+1
                                        endif
                                    endif
                                enddo
                            enddo
                        else
                            do l=1,i1-1
                                do m=1,lp
                                    if ((label(l,m).ne.c1(i,j)).and.(label(l,m).ne.c2(i,j)))then
                                        if(label(l,m).gt.c1(i,j))then
                                            label(l,m)=label(l,m)-1
                                        endif
                                    endif
                                enddo
                            enddo
                            do l=i1,i1
                                do m=1,j1-1
                                    if ((label(l,m).ne.c1(i,j)) .and.(label(l,m).ne.c2(i,j)))then
                                        if(label(l,m).gt.c1(i,j))then
                                            label(l,m)=label(l,m)-1
                                        endif
                                    endif
                                enddo
                            enddo
                        endif
                        
                    endif
                endif

                if((ab3.gt.0).and.(ab4.gt.0).and.(ab3.ne.ab4))then
                    l2=min(ab3,ab4)
                    if(latvalue(i,j).eq.0)then
                        label(i,j)=-l2
                    else
                        label(i,j)=l2
                    endif
                    i1=i
                    j1=j
                    if (l2.eq.ab3)then
                        do l=1,i1-1
                            do m=1,lp
                                if (label(l,m).eq.c4(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l2
                                    else
                                        label(l,m)=l2
                                    endif
                                endif
                            enddo
                        enddo
                        do l=i1,i1
                            do m=1,j1
                                if (label(l,m).eq.c4(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l2
                                    else
                                        label(l,m)=l2
                                    endif
                                endif
                            enddo
                        enddo
                    else
                        do l=1,i1
                            do m=1,j1
                                if(label(l,m).eq.c3(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l2
                                    else
                                        label(l,m)=l2
                                    endif
                                endif
                            enddo
                        enddo
                    endif
                endif

                if((ab1.gt.0).and.(ab3.gt.0).and.(ab1.ne.ab3))then
                    l3=min(ab1,ab3)
                    if(latvalue(i,j).eq.0)then
                        label(i,j)=-l3
                    else
                        label(i,j)=l3
                    endif
                    i1=i
                    j1=j
                    if (l1.eq.ab1)then
                        do l=1,i1-1
                            do m=1,lp
                                if (label(l,m).eq.c3(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l3
                                    else
                                        label(l,m)=l3
                                    endif
                                endif
                            enddo
                        enddo
                        do l=i1,i1
                            do m=1,j1
                                if (label(l,m).eq.c3(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l3
                                    else
                                        label(l,m)=l3
                                    endif
                                endif
                            enddo
                        enddo
                    else
                        do l=1,i1-1
                            do m=1,lp
                                if(label(l,m).eq.c1(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l3
                                    else
                                        label(l,m)=l3
                                    endif
                                endif
                            enddo
                        enddo
                        do l=i1,i1
                            do m=1,j1-1
                                if(label(l,m).eq.c1(i,j))then
                                    if(latvalue(i,j).eq.0)then
                                        label(l,m)=-l3
                                    else
                                        label(l,m)=l3
                                    endif
                                endif
                            enddo
                        enddo
                    endif
                endif

                if((c1(i,j).eq.0).and.(c2(i,j).eq.0).and.(c3(i,j).eq.0))then
                    if(c4(i,j).eq.0) then
                        if (latvalue(i,j).eq.0)then
                            i1=i
                            j1=j
                            m1= voidlabel
                            do l=1,i1-1
                                do m=1,lp
                                    if(label(l,m).lt.m1)then
                                        m1=label(l,m)
                                    endif
                                enddo
                            enddo
                            l=i1
                            do m=1,(j1-1)
                                if(label(l,m).lt.m1)then
                                    m1=label(l,m)
                                endif
                            enddo
                            voidlabel=m1
                            label(i,j)=voidlabel-1
                        else
                            i1=i
                            j1=j
                            m2= blacklabel
                            do l=1,i1-1
                                do m=1,lp
                                    if(label(l,m).gt.m2)then
                                        m2=label(l,m)
                                    endif
                                enddo
                            enddo
                            l=i1
                            do m=1,(j1-1)
                                if(label(l,m).gt.m2)then
                                    m2=label(l,m)
                                endif
                            enddo
                            blacklabel=m2
                            label(i,j)=blacklabel+1
                        endif
                    endif
                endif
            enddo
        enddo

        nb=0
        wl=0
        do i=1,lp
            do j=1,lp
                if(label(i,j).gt.nb)then
                    nb=nb+1
                endif
                if (label(i,j).lt.wl)then
                    wl=wl-1
                endif
            enddo
        enddo

        nw = abs(wl)
        nbt=nbt+nb
        nwt=nwt+nw

        do k=1,nb
            sizeb(k)=0
            do i=1,lp
                do j=1,lp
                    if(label(i,j).eq.k)then
                        sizeb(k)= sizeb(k)+1
                    endif
                enddo
            enddo
        enddo

        do n=1,nw
            sizew(n)=0
            do i=1,lp
                do j=1,lp
                    if(label(i,j).eq.(-n))then
                        sizew(n)=sizew(n)+1
                    endif
                enddo
            enddo
        enddo

        do k=1,nb
            clustersize(tot)=sizeb(k)
            trialno(tot)=itrial
            tot=tot+1
        enddo

    enddo  !trial ends

    nbf=nbt/trial
    nwf=nwt/trial
    eun=nbf-nwf
    eun1=real(eun)/real(lp*lp)
    
    chi = eun
    
    write(*,*)'Nb, Nw, eular no, normalized eular no'
    write(*,*)nbf,nwf,eun, eun1  ! Nb, Nw, eular no, normalized

    !do isz=1,lp*lp
    !    i=1
    !    tr(isz)=1
    !    do k =1,tot-1
    !        if (clustersize(k).eq.isz)then
    !            num(isz)=num(isz)+1
    !            trno(i)=trialno(k)
    !            i=i+1
    !            if((i.gt.2).and.(trno(i-1).ne.trno(i-2)))then
    !                tr(isz)=tr(isz)+1
    !            endif
    !        endif
    !    enddo
    !    num1(isz)=num(isz)/tr(isz)
    !enddo

end subroutine sq_chi

