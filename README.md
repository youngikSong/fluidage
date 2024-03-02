# fluidage
fluid age code for Openfoam

timeTransport 

fluid-age transport equation. Uses variable s, which is actually same from original scalartransportFoam, so watch out if you add another

mfractionTransport

mass fraction transport equation, uses variable zeta

avgageTransport

mass-weighted stream-age equation, uses variable cphi as acronym for capital phi, also uses variable zeta so should be used together with mfraction