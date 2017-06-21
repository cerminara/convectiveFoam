How to use convectiveFoam?

- clone the source files in your preferred directory:
    git clone https://github.com/cerminara/convectiveFoam.git

- Be sure the environment variables of OpenFoam 4.x have been
  loaded in your terminal window

- enter in the convectiveFoam folder and install it launching 
  the script:
    ./Allwmake

- try it with the singlePlume case, in your tutorial folder, by
  using the following commands:
    cd $WM_PROJECT_USER_DIR/tutorials/convectiveFoam/singlePlume/Pr1_Ra1e6_limLin
    ./Allrun
