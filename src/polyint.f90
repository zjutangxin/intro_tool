! test Polynomial Int

program polytest

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: nb = 20
    real(dp), parameter :: bmin = 0.0_dp, bmax = 0.8_dp

    integer :: indb
    real(dp), dimension(nb) :: bvec
    real(dp), dimension(nb,nb) :: DebPol1, DebPol2
    real(dp), dimension(nb,nb) :: DebPol1fit, DebPol2fit
    real(dp), dimension(7) :: coevec1, coevec2

    do indb = 1, nb
        bvec(indb) = bmin+(bmax-bmin)*(real(indb-1)/real(nb-1))
    end do
    
    open(1,file='./DebPol1_50.txt',form='formatted')
    do indb = 1,nb
        read(1,'(20ES14.6)') DebPol1(indb,:)
    end do
    close(1)

    open(1,file='./DebPol2_50.txt',form='formatted')
    do indb = 1,nb
        read(1,'(20ES14.6)') DebPol2(indb,:)
    end do
    close(1)

    call PolynomialInt(bVec,Nb,bVec,Nb,DebPol1,DebPol1fit,CoeVec1)
    call PolynomialInt(bVec,Nb,bVec,Nb,DebPol2,DebPol2fit,CoeVec2)

    write (*,'(A10,7f10.6)') 'coev1 = ', coevec1
    write (*,'(A10,7f10.6)') 'coev2 = ', coevec2

    open(1,file='./DebPol1fit.txt',form='formatted')
    do indb = 1,nb
        write(1,'(20ES14.6)') DebPol1fit(indb,:)
    end do
    close(1)

    open(1,file='./DebPol2fit.txt',form='formatted')
    do indb = 1,nb
        write(1,'(20ES14.6)') DebPol2fit(indb,:)
    end do
    close(1)

end program polytest

subroutine PolynomialInt(xVec,nxR,yVec,nyR,zMx,zzMx,CoeVec)
    ! polynomial fit
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    integer, intent(in) :: nxr,nyr
    real(dp), dimension(nxr), intent(in) :: xvec
    real(dp), dimension(nyr), intent(in) :: yvec
    real(dp), dimension(nxr,nyr), intent(in) :: zmx
    real(dp), dimension(nxr,nyr), intent(out) :: zzmx
    real(dp), dimension(7), intent(out) :: coevec

    integer :: indx,indy,indxy,info
    integer, parameter :: block_size = 32
    integer, parameter :: lwork = 7+7*block_size
    real(dp), dimension(nxr*nyr) :: DepVarVec
    real(dp), dimension(nxr*nyr,7) :: IndVarMx, IndVarMx_pre
    real(dp), dimension(lwork) :: work
        
    do indx = 1, nxR
        do indy = 1, nyR
            indxy=(indx-1)*nyR+indy
            IndVarMx(indxy,1) = 1.0_dp
            IndVarMx(indxy,2) = xVec(indx)
            IndVarMx(indxy,3) = yVec(indy)
            IndVarMx(indxy,4) = xVec(indx)*yVec(indy)
            IndVarMx(indxy,5) = xVec(indx)**2.0_dp
            IndVarMx(indxy,6) = yVec(indy)**2.0_dp
            IndVarMx(indxy,7) = xVec(indx)**2.0_dp*yVec(indy)**2.0_dp
            DepVarVec(indxy) = zMx(indx,indy)
        end do
    end do
    IndVarMx_pre = IndVarMx
        ! syntax
        ! DGELS(TRANS,M,N,NRHS,A,LDA,B,LDB,WORK,LWORK,INFO)
    call DGELS('N',nxr*nyr,7,1,IndVarMx,nxr*nyr,DepVarVec,nxr*nyr,work,lwork,info)
    coevec = DepVarVec
        ! CALL DRLSE(nxR*nyR,DepVarVec,6,IndVarMx(:,2:7),nxR*nyR,1,coeVec,SST,SSE)
    do indx = 1, nxR
        do indy = 1, nyR
            indxy = (indx-1)*nyR+indy
            zzMx(indx,indy) = dot_product(IndVarMx_pre(indxy,:),coevec)
        end do
    end do
end subroutine PolynomialInt