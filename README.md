# flex-fea

This is a tool for designing flexible structures using finite element modeling.

# requirements

You will need a relatively new version of matlab (>=R2017b) because the syntax for setting up PDEs has changed (see the [PDE toolbox changelog](https://www.mathworks.com/help/pde/release-notes.html)).

# demo code

A sample model can be created and shown with the following commands:
```
stl_filename = 'ultra-amp-4mm.stl'
flex = UltraFlexure(stl_filename);
flex.solveModal();
flex.showModes();
modes_to_show = 4:9;
flex.generateFrames(modes_to_show);
flex.showMovie();
```
