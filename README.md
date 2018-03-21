# flex-fea

This is a tool for designing flexible structures using finite element modeling.

# requirements

You will need a relatively new version of matlab (>=R2017b) because the syntax for setting up PDEs has changed (see the [PDE toolbox changelog](https://www.mathworks.com/help/pde/release-notes.html)).

# demo code

A sample model can be created and shown with the following commands:
```
flex = ParallelFlex(1,15,45,15,40,0)
flex.showDeflectionY()
```
