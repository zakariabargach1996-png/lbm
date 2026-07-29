module milestone4
  ! Minimal library entry point retained from the fpm project scaffold.  The
  ! actual simulation executable is implemented in app/main.f90.
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    ! Print the library scaffold greeting.
    print *, "Hello, milestone4!"
  end subroutine say_hello
end module milestone4
