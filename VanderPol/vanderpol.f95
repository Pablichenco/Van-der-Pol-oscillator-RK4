Program vanderpol

    implicit none
    real, dimension(2) :: x
    real :: t
    real :: h
    integer :: n
    integer :: nf ! n final
    character (len=100):: archv!el nombre del archivo
    !inicializar variables
    n=2
    nf=5000
    t=0.0
    h=0.1

    !pregunta las condiciones iniciamen manteniendo a to=0 fijo
    print*, 'Condiciones iniciales x0 , dx/dt0 (v0), recuerde que t_0=0.0:'
    read*,x

    !pregunta el nombre del archivo
    print*,'Escribir el nombre del archivo CON la extencion corespondiente,  (ejemplo: datos.txt)'
    read*,archv


    open (unit=4, File =archv, status ='unknown')

    !llamar a la subrutina de Runge Kutta 4
    call RungeKutta4(t,x,h,n,nf)
    stop

    close(4)

end Program vanderpol

subroutine xpol(x,k)
    real :: M ! M=µ
    real, dimension(2) :: x,k
    M=2.5         ! valor de µ para modificar
    k(1)=x(2)
    k(2)=M*(1-x(1)**2)*x(2)-x(1)
    return
end subroutine xpol

!subrutina de Runge-Kutta
subroutine RungeKutta4(t,x,h,n,nf)!n,nstep=np,k=l,I=j

    real ::h,h2,t,ti
    integer :: n,nf,i,j
    real, dimension(10) :: x,y,k1,k2,k3,k4

    !iniciar h2=0.5h ti=t
    h2=0.5*h !h2=h/2
    ti=t

    write(*,*) t,x(1),x(2)
    write(4,*) t,x(1),x(2)

    do  i=1,nf !primer do #1

        call xpol(x,k1) !inicia k1

        do  j=1,n !segundo do #2
            y(j)=x(j)+h2*k1(j)
        end do !fin do#2

         call xpol(y,k2)

        do  j=1,n  !tercer do #3
            y(j)=x(j)+h2*k2(j) !inicia k2
        end do!fin do#3

        call xpol(y,k3)

        do  j=1,n !cuarto do #4
            y(j)=x(j)+h2*k3(j) !inicial k3
        end do!fin do#4

        call xpol(y,k4)

        do  j=1,n !quinto do #5
            x(j)=x(j)+h*( k1(j) + 2.0*( k2(j) + k3(j) ) + k4(j) )/6.0
        end do!fin do#5

        t=ti+h*real(i)
        write(*,*)t,(x(j),j=1,n)
        write(4,*) t,(x(j),j=1,n)
        write(4,*)t,x(2)
        end do!fin do#1
        return
    close(4)


end subroutine RungeKutta4

