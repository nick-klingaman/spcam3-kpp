	  �3  d   k820309              15.0        ޽W                                                                                                           
       /home2/n02/n02/pappas/src/cam3_sp/models/lnd/clm2/src/main/lp_coupling.F90 LP_COUPLING              LP_COUPLING_INIT LP_COUPLING_FINALIZE ALLTOALL_CLUMP_TO_CHUNK_INIT ALLTOALL_CLUMP_TO_CHUNK ALLTOALL_CHUNK_TO_CLUMP                                                     
       MPIR8 MPICOM          @       �   @                              
       NPES                                                     
       IAM                      @                              
       GET_NCLUMPS GET_CLUMP_OWNER_ID GET_CLUMP_NCELLS_ID GET_CLUMP_COORD_ID GET_CLUMP_GCELL_INFO                      @                              
       GET_CHUNK_COORD_OWNER_P                  � @                              
       R8 SHR_KIND_R8            @@                                                      @@                                                                                      	            %         @                               
                            %         @                                                          #CID              
                                             %         @                                                          #CID              
                                             #         @                                                      #CID    #NCELLS    #LONS    #LATS              
                                                       
                                                                                                                p          5 O p            5 O p                                                                                         	    p          5 O p            5 O p                          #         @                                                      #CID    #CELL    #GI              
                                                       
                                                                                                                                                          #         @                                                       #LP_COUPLING_INIT%SUM    #LP_COUPLING_INIT%MAXVAL                  @                                SUM               @                                MAXVAL #         @                                                       #         @                                                    #ALLTOALL_CLUMP_TO_CHUNK_INIT%NPES    #ALLTOALL_CLUMP_TO_CHUNK_INIT%ENDCHUNK    #ALLTOALL_CLUMP_TO_CHUNK_INIT%BEGCHUNK     #ALLTOALL_CLUMP_TO_CHUNK_INIT%ENDCHUNK !   #ALLTOALL_CLUMP_TO_CHUNK_INIT%BEGCHUNK "   #ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM #   #ALLTOALL_CLUMP_TO_CHUNK_INIT%SQRT 1   #SRFFLX2D 2                                                          @                           #     '�                   #NAME $   #ASDIR %   #ASDIF &   #ALDIR '   #ALDIF (   #LWUP )   #LHF *   #SHF +   #WSX ,   #WSY -   #TREF .   #TS /   #CFLX 0                �                              $                                        �                              %                             
  p          p            p                                       �                              &            X                 
  p          p            p                                       �                              '            �                 
  p          p            p                                       �                              (            �                 
  p          p            p                                       �                              )                            
  p          p            p                                       �                              *            X                
  p          p            p                                       �                              +            �                
  p          p            p                                       �                              ,            �             	   
  p          p            p                                       �                              -                         
   
  p          p            p                                       �                              .            X                
  p          p            p                                       �                              /            �                
  p          p            p                                       �                              0            �                
  p 	         p          p            p          p                                                                                                                                                                                                                                     !                                                       "                          @                           1     SQRT          
D     �                           2             �            5 r "     & 5 r "   5 r !         5 r !   5 r "   p                          #ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM #   #         @                                   3                  #ALLTOALL_CLUMP_TO_CHUNK%NPES 4   #ALLTOALL_CLUMP_TO_CHUNK%ENDCHUNK 5   #ALLTOALL_CLUMP_TO_CHUNK%BEGCHUNK 6   #ALLTOALL_CLUMP_TO_CHUNK%ENDCHUNK 7   #ALLTOALL_CLUMP_TO_CHUNK%BEGCHUNK 8   #ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM 9   #SRFFLX2D G                                                                                                                           @                           9     '�                   #NAME :   #ASDIR ;   #ASDIF <   #ALDIR =   #ALDIF >   #LWUP ?   #LHF @   #SHF A   #WSX B   #WSY C   #TREF D   #TS E   #CFLX F                �                              :                                        �                              ;                             
  p          p            p                                       �                              <            X                 
  p          p            p                                       �                              =            �                 
  p          p            p                                       �                              >            �                 
  p          p            p                                       �                              ?                            
  p          p            p                                       �                              @            X                
  p          p            p                                       �                              A            �                
  p          p            p                                       �                              B            �             	   
  p          p            p                                       �                              C                         
   
  p          p            p                                       �                              D            X                
  p          p            p                                       �                              E            �                
  p          p            p                                       �                              F            �                
  p 	         p          p            p          p                                                                   4                                                     5                                                     6                                                       7                                                       8                     
D     �                           G             �            5 r 8     & 5 r 8   5 r 7         5 r 7   5 r 8   p                          #ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM 9   #         @                                   H                 #ALLTOALL_CHUNK_TO_CLUMP%NPES I   #ALLTOALL_CHUNK_TO_CLUMP%ENDCHUNK J   #ALLTOALL_CHUNK_TO_CLUMP%BEGCHUNK K   #ALLTOALL_CHUNK_TO_CLUMP%ENDCHUNK L   #ALLTOALL_CHUNK_TO_CLUMP%BEGCHUNK M   #ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE N   #ALLTOALL_CHUNK_TO_CLUMP%SQRT a   #SRF_STATE b                                                     @                           N     '@                   #TBOT O   #ZBOT P   #UBOT Q   #VBOT R   #QBOT S   #PBOT T   #FLWDS U   #PRECSC V   #PRECSL W   #PRECC X   #PRECL Y   #SOLL Z   #SOLS [   #SOLLD \   #SOLSD ]   #SRFRAD ^   #THBOT _   #TSSUB `                �                              O                              
  p          p            p                                       �                              P            @                 
  p          p            p                                       �                              Q            �                 
  p          p            p                                       �                              R            �                 
  p          p            p                                       �                              S                             
  p          p            p                                       �                              T            @                
  p          p            p                                       �                              U            �                
  p          p            p                                       �                              V            �                
  p          p            p                                       �                              W                          	   
  p          p            p                                       �                              X            @             
   
  p          p            p                                       �                              Y            �                
  p          p            p                                       �                              Z            �                
  p          p            p                                       �                              [                             
  p          p            p                                       �                              \            @                
  p          p            p                                       �                              ]            �                
  p          p            p                                       �                              ^            �                
  p          p            p                                       �                              _                             
  p          p            p                                       �                              `             @                
  p 	         p          p            p          p                                                                   I                                                     J                                                     K                                                       L                                                       M                          @                           a     SQRT          
      �                           b             @           5 r M     & 5 r M   5 r L         5 r L   5 r M   p                          #ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE N      �   _      fn#fn !   �   �   b   uapp(LP_COUPLING    �  M   J  MPISHORTHAND    �  E   J  SPMD_DYN      D   J  PMGRID    X  �   J  LND_GRID    �  X   J  PHYS_GRID    K  O   J  SHR_KIND_MOD #   �  @       MPIR8+MPISHORTHAND $   �  @       MPICOM+MPISHORTHAND      @       NPES+SPMD_DYN %   Z  P       GET_NCLUMPS+LND_GRID ,   �  Y       GET_CLUMP_OWNER_ID+LND_GRID 0     @   a   GET_CLUMP_OWNER_ID%CID+LND_GRID -   C  Y       GET_CLUMP_NCELLS_ID+LND_GRID 1   �  @   a   GET_CLUMP_NCELLS_ID%CID+LND_GRID ,   �  q       GET_CLUMP_COORD_ID+LND_GRID 0   M  @   a   GET_CLUMP_COORD_ID%CID+LND_GRID 3   �  @   a   GET_CLUMP_COORD_ID%NCELLS+LND_GRID 1   �  �   a   GET_CLUMP_COORD_ID%LONS+LND_GRID 1   q  �   a   GET_CLUMP_COORD_ID%LATS+LND_GRID .     c       GET_CLUMP_GCELL_INFO+LND_GRID 2   x  @   a   GET_CLUMP_GCELL_INFO%CID+LND_GRID 3   �  @   a   GET_CLUMP_GCELL_INFO%CELL+LND_GRID 1   �  @   a   GET_CLUMP_GCELL_INFO%GI+LND_GRID    8	  @       NPES+SPMD_DYN !   x	         LP_COUPLING_INIT %   �	  <      LP_COUPLING_INIT%SUM (   3
  ?      LP_COUPLING_INIT%MAXVAL %   r
  H       LP_COUPLING_FINALIZE -   �
  �      ALLTOALL_CLUMP_TO_CHUNK_INIT @   ]  �      ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM+COMSRF E   -  P   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%NAME+COMSRF F   }  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%ASDIR+COMSRF F     �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%ASDIF+COMSRF F   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%ALDIR+COMSRF F   Q  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%ALDIF+COMSRF E   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%LWUP+COMSRF D   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%LHF+COMSRF D   %  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%SHF+COMSRF D   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%WSX+COMSRF D   ]  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%WSY+COMSRF E   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%TREF+COMSRF C   �  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%TS+COMSRF E   1  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX_PARM%CFLX+COMSRF @   �  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%NPES+SPMD_DYN=NPES F   -  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%ENDCHUNK+PPGRID=ENDCHUNK F   m  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%BEGCHUNK+PPGRID=BEGCHUNK =   �  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%ENDCHUNK+PPGRID =   �  @     ALLTOALL_CLUMP_TO_CHUNK_INIT%BEGCHUNK+PPGRID 2   -  =      ALLTOALL_CLUMP_TO_CHUNK_INIT%SQRT 6   j  �   a   ALLTOALL_CLUMP_TO_CHUNK_INIT%SRFFLX2D (   \  �      ALLTOALL_CLUMP_TO_CHUNK ;   �  �      ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM+COMSRF @   �  P   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%NAME+COMSRF A     �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%ASDIR+COMSRF A   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%ASDIF+COMSRF A   S  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%ALDIR+COMSRF A   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%ALDIF+COMSRF @   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%LWUP+COMSRF ?   '  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%LHF+COMSRF ?   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%SHF+COMSRF ?   _  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%WSX+COMSRF ?   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%WSY+COMSRF @   �  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%TREF+COMSRF >   3   �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%TS+COMSRF @   �   �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX_PARM%CFLX+COMSRF ;   �!  @     ALLTOALL_CLUMP_TO_CHUNK%NPES+SPMD_DYN=NPES A   �!  @     ALLTOALL_CLUMP_TO_CHUNK%ENDCHUNK+PPGRID=ENDCHUNK A   "  @     ALLTOALL_CLUMP_TO_CHUNK%BEGCHUNK+PPGRID=BEGCHUNK 8   K"  @     ALLTOALL_CLUMP_TO_CHUNK%ENDCHUNK+PPGRID 8   �"  @     ALLTOALL_CLUMP_TO_CHUNK%BEGCHUNK+PPGRID 1   �"  �   a   ALLTOALL_CLUMP_TO_CHUNK%SRFFLX2D (   �#  ~      ALLTOALL_CHUNK_TO_CLUMP =   6%       ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE+COMSRF B   G&  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%TBOT+COMSRF B   �&  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%ZBOT+COMSRF B   '  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%UBOT+COMSRF B   (  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%VBOT+COMSRF B   �(  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%QBOT+COMSRF B   S)  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PBOT+COMSRF C   �)  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%FLWDS+COMSRF D   �*  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PRECSC+COMSRF D   '+  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PRECSL+COMSRF C   �+  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PRECC+COMSRF C   _,  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%PRECL+COMSRF B   �,  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SOLL+COMSRF B   �-  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SOLS+COMSRF C   3.  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SOLLD+COMSRF C   �.  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SOLSD+COMSRF D   k/  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%SRFRAD+COMSRF C   0  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%THBOT+COMSRF C   �0  �   a   ALLTOALL_CHUNK_TO_CLUMP%SURFACE_STATE%TSSUB+COMSRF ;   _1  @     ALLTOALL_CHUNK_TO_CLUMP%NPES+SPMD_DYN=NPES A   �1  @     ALLTOALL_CHUNK_TO_CLUMP%ENDCHUNK+PPGRID=ENDCHUNK A   �1  @     ALLTOALL_CHUNK_TO_CLUMP%BEGCHUNK+PPGRID=BEGCHUNK 8   2  @     ALLTOALL_CHUNK_TO_CLUMP%ENDCHUNK+PPGRID 8   _2  @     ALLTOALL_CHUNK_TO_CLUMP%BEGCHUNK+PPGRID -   �2  =      ALLTOALL_CHUNK_TO_CLUMP%SQRT 2   �2  �   a   ALLTOALL_CHUNK_TO_CLUMP%SRF_STATE 