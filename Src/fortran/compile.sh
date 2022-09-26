

  gfortran -o protein.exe ddot.f90 dpmeps.f90  imported.f90 timer.f90 protein.f90 -Og  -fcheck=all -fbacktrace
  gfortran -o protein.exe ddot.f90 dpmeps.f90  imported.f90 timer.f90 protein.f90


  cp protein.exe ../testing

