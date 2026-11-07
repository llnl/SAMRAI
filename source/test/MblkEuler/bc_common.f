c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   Define some boundary condition constants.
c

c
c     define faces, edges and corners that match up to 
c     apputils/boundary/CartesianBoundaryDefines.h
c

        integer XLO      
        integer XHI      
        integer YLO      
        integer YHI      
        integer ZLO      
        integer ZHI      

        integer Y0Z0     
        integer Y1Z0     
        integer Y0Z1     
        integer Y1Z1     
        integer X0Z0     
        integer X0Z1     
        integer X1Z0     
        integer X1Z1     
        integer X0Y0     
        integer X1Y0     
        integer X0Y1     
        integer X1Y1     

        integer X0Y0Z0   
        integer X1Y0Z0   
        integer X0Y1Z0   
        integer X1Y1Z0   
        integer X0Y0Z1   
        integer X1Y0Z1   
        integer X0Y1Z1   
        integer X1Y1Z1   

        parameter( XLO = 0)
        parameter( XHI = 1)
        parameter( YLO = 2)
        parameter( YHI = 3)
        parameter( ZLO = 4)
        parameter( ZHI = 5)

        parameter( Y0Z0 = 0)  
        parameter( Y1Z0 = 1)
        parameter( Y0Z1 = 2)
        parameter( Y1Z1 = 3)
        parameter( X0Z0 = 4)
        parameter( X0Z1 = 5)
        parameter( X1Z0 = 6)
        parameter( X1Z1 = 7)
        parameter( X0Y0 = 8)
        parameter( X1Y0 = 9)
        parameter( X0Y1 = 10)
        parameter( X1Y1 = 11)

        parameter( X0Y0Z0 = 0)
        parameter( X1Y0Z0 = 1)
        parameter( X0Y1Z0 = 2)
        parameter( X1Y1Z0 = 3)
        parameter( X0Y0Z1 = 4)
        parameter( X1Y0Z1 = 5)
        parameter( X0Y1Z1 = 6)
        parameter( X1Y1Z1 = 7)
c
c     end of include file
c
