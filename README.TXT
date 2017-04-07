                          

          P A F : A software tool to estimate free-geometry extended bodies of anomalous pressure from surface deformation data.


          by  A.G. Camacho (1), J. Fernández (1*), F. Cannavó (2)
          (1) Institute of Geosciences (CSIC-UCM), Madrid, Spain
          (2) Osservatorio Etneo, Istituto Nazionale di Geofisica e Vulcanologia, Catania, Italy. 
           *Corresponding author:e-mail: jft@mat.ucm.es; phone: +34-913944632


The software package PAF carries out inversions of surface deformation data (InSAR, GPS, and terrestrial data) as produced by 3D free-geometry extended bodies with anomalous pressure changes. The anomalous structures are described as aggregation of elementary cells (point sources) in an elastic half space. This software can be used to invert large data sets (e.g. InSAR data sets) or with data coming from small control networks (e.g. GPS monitoring data), mainly in volcanic areas to estimate the expected causative pressure bodies.


This package includes:


    (a) A User Manual: PAFmanual.txt


    (b) Two executable files (if files appears as  *-exe, please rename as *.exe to execute): 

                  ConfigPAF.exe
 
                  InverPAF.exe 


    (c) The corresponding codes for the previous executable files:

                  ConfigPAF.for
                  
                  InverPAF.exe                                


    (c) Files for an application example (*)


                  Data files:   DeforData.txt , and map.bln  (the second one,only for pictures)
 
                  Intermediate files:   CellsConfig.txt 

                  Results files:    ModPAF.exe 


(*) This application case correspond to a simulated deep vertical ellipsoid for negative pressure, located at a mean depth 3 km b.s.l.
We have taken the following parameters:
          Depth to model bottom: 6000 m
          Mean side of cells:     150 m
          Smoothing coefficient:    2
          Press contrast:          60  MPa
         Significance limit:       3

