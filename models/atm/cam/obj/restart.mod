	  �4  �   k820309              15.0        ٽW                                                                                                           
       /home2/n02/n02/pappas/src/cam3_sp_intel/models/atm/cam/src/control/restart.F90 RESTART              WRITE_RESTART READ_RESTART SET_RESTART_FILEPATH                                                     
       MASTERPROC PLEV PLEVP PLOND PLAT                                                     
       NLON WNUMMAX                      @                              
       PUTFIL GETFIL OPNFIL                      @                              
       GET_ARCHIVEDIR          @       �                                  
       MPICOM MPIR8 MPIINT MPILOG                                                     
       R8 SHR_KIND_R8               � @                              '                     #DAY    #TOD 	                � D                                                            � D                             	                          #ESMF_TOD 
                  � @                         
     '                    #TYPE    #SEC    #MSEC                 � D                                                            � D                                                           � D                                                            � @                              '                    #CALENDAR    #YEAR    #MONTH    #DAY    #TOD    #JULIANDAY    #DAYOFYEAR                 � D                                  �                      #ESMF_CALENDAR                   � @                              '�                    #TYPE    #DIM    #DIMRUNNINGSUM    #DIY                 � D                                                            � D                                                           p          p            p                                       � D                                        p                   p          p            p                                       � D                                 �                          � D                                 �                          � D                                 �                          � D                                 �                          � D                                         �              #ESMF_TOD 
                � D                                                          � D                                                           � @                               '�                   #NSTEP    #STEPSIZE    #STARTDATE    #STOPDATE    #BASEDATE     #CURRDATE !   #PREVDATE "                � D                                                            � D                                                        #ESMF_TIME                 � D                                         (              #ESMF_DATE                 � D                                         H             #ESMF_DATE                 � D                                          h             #ESMF_DATE                 � D                             !            �             #ESMF_DATE                 � D                             "            �             #ESMF_DATE    #         @                                  #                   #PUTFIL%TRIM $   #LOCFN %   #MSSFPN &   #PASS '   #IRT (   #LREMOV )                 @                           $     TRIM           
                                 %                    1           
                                 &                    1           
                                 '                    1           
                                  (                     
                                  )           #         @                                  *                   #GETFIL%LEN_TRIM +   #GETFIL%TRIM ,   #GETFIL%PRESENT -   #FULPATH .   #LOCFN /   #IFLAG 0                 @                           +     LEN_TRIM               @                           ,     TRIM               @                           -     PRESENT           
                                 .                    1                                           /                     1           
                                 0           #         @                                  1                   #OPNFIL%LEN_TRIM 2   #OPNFIL%TRIM 3   #OPNFIL%PRESENT 4   #LOCFN 5   #IUN 6   #FORM 7   #STATUS 8                 @                           2     LEN_TRIM               @                           3     TRIM               @                           4     PRESENT           
                                 5                    1           
                                  6                     
                                 7                                     
                                8                    1 $         @                                 9                          #GET_ARCHIVEDIR%TRIM :   #TYPE ;                         @                           :     TRIM           
                                 ;                    1            @@                               <                       @@                               =                       @@                               >                       @@                               ?            #         @                                   @                   #WRITE_RESTART%NUMLATS A   #WRITE_RESTART%NPES B   #WRITE_RESTART%TRIM C                                             A                                                     B                          @                           C     TRIM #         @                                   D     	               #READ_RESTART%CHUNK_BUF_NRECS E   #READ_RESTART%BLOCK_BUF_NRECS F   #READ_RESTART%NLCOLS G   #READ_RESTART%NGCOLS H   #READ_RESTART%NCHUNKS I   #READ_RESTART%ENDCHUNK J   #READ_RESTART%BEGCHUNK K   #READ_RESTART%NUMLATS L   #READ_RESTART%NPES M                                            E                                                     F                                                     G                                                     H                                                     I                                                     J                                                     K                                                      L                                                     M            #         @                                   N                   #SET_RESTART_FILEPATH%LEN_TRIM O   #SET_RESTART_FILEPATH%TRIM P   #RGPATH Q                 @                           O     LEN_TRIM               @                           P     TRIM           
  @                              Q                    1                @                           R     �       "              #ITSST S   #NSREST T   #IRADSW U   #IRADLW V   #IRADAE W   #NREFRQ X   #ANNCYC Y   #NLEND Z   #NLRES [   #NLHST \   #LBRNCH ]   #AERES ^   #OZNCYC _   #SSTCYC `   #ICECYC a   #ADIABATIC b   #FLXAVE c   #IDEAL_PHYS d   #NSPLIT e   #IORD f   #JORD g   #KORD h   #USE_ETA i   #AQUA_PLANET j   #DORAMP_SO4 k   #DORAMP_SCON l   #FULLGRID m   #PRINT_STEP_COST n   #DOABSEMS o   #DOSW p   #DOLW q   #INDIRECT r   #SOM_CONSCHK_FRQ s   #ICE_CONSCHK_FRQ t             �   @        �                   S                              �   @        �                   T                             �   @        �                   U                             �   @        �                   V                             �   @        �                   W                             �   @        �                   X                             �   @        �                   Y                             �   @        �                   Z                             �   @        �                   [                              �   @        �                   \     $                        �   @        �                   ]     (                        �@ @@        �                   ^     ,                        �   @        �                   _     0                        �   @        �                   `     4                        �   @        �                   a     8                        �   @        �                   b     <                        �   @        �                   c     @                        �   @        �                   d     D                        �   @        �                   e     H                        �   @        �                   f     L                        �   @        �                   g     P                        �   @        �                   h     T                        �   @        �                   i     X                        �   @        �                   j     \                        �   @        �                   k     `                        �   @        �                   l     d                        �   @        �                   m     h                        �   @        �                   n     l                        �   @        �                   o     p                        �   @        �                   p     t                        �   @        �                   q     x                        �   @        �                   r     |                        �   @        �                   s     �                        �   @        �                   t     �                             @                           u                          #DIVDAMPN v   #PRECC_THRESH w   #PRECL_THRESH x             �   @       �                   v             
                 �   @       �                   w            
                 �   @       �                   x            
                      @                           y                           #NSDS z   #NRG {   #NRG2 |   #NCID_INI }   #NCID_OZ ~   #NCID_SST    #NCID_TRC �   #LUHREST �             �   @        �                   z                              �@ @@        �                   {                             �@ @@        �                   |                             �   @        �                   }                             �   @        �                   ~                             �   @        �                                                �   @        �                   �                             �  @@        �                   �                                  @                           �     �                    #HYAI �   #HYAM �   #HYBI �   #HYBM �   #HYBD �   #HYPI �   #HYPM �   #HYPD �   #PS0 �   #PSR �   #PRSFAC �   #NPRLEV �             �@ @@       �                   �                           
      p          p            p                                            �@ @@       �                   �                   �       
      p          p            p                                            �@ @@       �                   �                   �      
      p          p            p                                            �@ @@       �                   �                   �      
      p          p            p                                            �   @       �                   �                   �      
      p          p            p                                            �   @       �                   �                   �      
      p          p            p                                            �   @       �                   �                   �      
      p          p            p                                            �   @       �                   �                   �      
      p          p            p                                            �   @       �                   �     �      
                 �   @       �                   �     �      
                 �   @       �                   �     �      
                 �   @        �                   �     �                            @                           �                          #EPS �             �@ @@       �                   �             
          �   _      fn#fn    �   @   b   uapp(RESTART    ?  a   J  PMGRID    �  M   J  RGRID    �  U   J  IOFILEMOD    B  O   J  FILENAMES    �  [   J  MPISHORTHAND    �  O   J  SHR_KIND_MOD '   ;  b       ESMF_TIME+ESMF_TIMEMOD /   �  H   %   ESMF_TIME%DAY+ESMF_TIMEMOD=DAY /   �  ^   %   ESMF_TIME%TOD+ESMF_TIMEMOD=TOD %   C  m      ESMF_TOD+ESMF_TODMOD /   �  H   %   ESMF_TOD%TYPE+ESMF_TODMOD=TYPE -   �  H   %   ESMF_TOD%SEC+ESMF_TODMOD=SEC /   @  H   %   ESMF_TOD%MSEC+ESMF_TODMOD=MSEC '   �  �       ESMF_DATE+ESMF_DATEMOD 9   +  c   %   ESMF_DATE%CALENDAR+ESMF_DATEMOD=CALENDAR /   �        ESMF_CALENDAR+ESMF_CALENDARMOD 9     H   %   ESMF_CALENDAR%TYPE+ESMF_CALENDARMOD=TYPE 7   U  �   %   ESMF_CALENDAR%DIM+ESMF_CALENDARMOD=DIM K   �  �   %   ESMF_CALENDAR%DIMRUNNINGSUM+ESMF_CALENDARMOD=DIMRUNNINGSUM 7   �  H   %   ESMF_CALENDAR%DIY+ESMF_CALENDARMOD=DIY 1   �  H   %   ESMF_DATE%YEAR+ESMF_DATEMOD=YEAR 3   	  H   %   ESMF_DATE%MONTH+ESMF_DATEMOD=MONTH /   e	  H   %   ESMF_DATE%DAY+ESMF_DATEMOD=DAY /   �	  ^   %   ESMF_DATE%TOD+ESMF_DATEMOD=TOD ;   
  H   %   ESMF_DATE%JULIANDAY+ESMF_DATEMOD=JULIANDAY ;   S
  H   %   ESMF_DATE%DAYOFYEAR+ESMF_DATEMOD=DAYOFYEAR -   �
  �       ESMF_TIMEMGR+ESMF_TIMEMGRMOD 9   K  H   %   ESMF_TIMEMGR%NSTEP+ESMF_TIMEMGRMOD=NSTEP ?   �  _   %   ESMF_TIMEMGR%STEPSIZE+ESMF_TIMEMGRMOD=STEPSIZE A   �  _   %   ESMF_TIMEMGR%STARTDATE+ESMF_TIMEMGRMOD=STARTDATE ?   Q  _   %   ESMF_TIMEMGR%STOPDATE+ESMF_TIMEMGRMOD=STOPDATE ?   �  _   %   ESMF_TIMEMGR%BASEDATE+ESMF_TIMEMGRMOD=BASEDATE ?     _   %   ESMF_TIMEMGR%CURRDATE+ESMF_TIMEMGRMOD=CURRDATE ?   n  _   %   ESMF_TIMEMGR%PREVDATE+ESMF_TIMEMGRMOD=PREVDATE !   �  �       PUTFIL+IOFILEMOD +   \  =      PUTFIL%TRIM+IOFILEMOD=TRIM '   �  L   a   PUTFIL%LOCFN+IOFILEMOD (   �  L   a   PUTFIL%MSSFPN+IOFILEMOD &   1  L   a   PUTFIL%PASS+IOFILEMOD %   }  @   a   PUTFIL%IRT+IOFILEMOD (   �  @   a   PUTFIL%LREMOV+IOFILEMOD !   �  �       GETFIL+IOFILEMOD 3   �  A      GETFIL%LEN_TRIM+IOFILEMOD=LEN_TRIM +   �  =      GETFIL%TRIM+IOFILEMOD=TRIM 1      @      GETFIL%PRESENT+IOFILEMOD=PRESENT )   `  L   a   GETFIL%FULPATH+IOFILEMOD '   �  L   a   GETFIL%LOCFN+IOFILEMOD '   �  @   a   GETFIL%IFLAG+IOFILEMOD !   8  �       OPNFIL+IOFILEMOD 3   �  A      OPNFIL%LEN_TRIM+IOFILEMOD=LEN_TRIM +   %  =      OPNFIL%TRIM+IOFILEMOD=TRIM 1   b  @      OPNFIL%PRESENT+IOFILEMOD=PRESENT '   �  L   a   OPNFIL%LOCFN+IOFILEMOD %   �  @   a   OPNFIL%IUN+IOFILEMOD &   .  P   a   OPNFIL%FORM+IOFILEMOD (   ~  L   a   OPNFIL%STATUS+IOFILEMOD )   �  {       GET_ARCHIVEDIR+FILENAMES 3   E  =      GET_ARCHIVEDIR%TRIM+FILENAMES=TRIM .   �  L   a   GET_ARCHIVEDIR%TYPE+FILENAMES $   �  @       MPICOM+MPISHORTHAND #     @       MPIR8+MPISHORTHAND $   N  @       MPIINT+MPISHORTHAND $   �  @       MPILOG+MPISHORTHAND    �  �       WRITE_RESTART 5   a  @     WRITE_RESTART%NUMLATS+PMGRID=NUMLATS 1   �  @     WRITE_RESTART%NPES+SPMD_DYN=NPES #   �  =      WRITE_RESTART%TRIM      ?      READ_RESTART G   ]  @     READ_RESTART%CHUNK_BUF_NRECS+PHYS_GRID=CHUNK_BUF_NRECS G   �  @     READ_RESTART%BLOCK_BUF_NRECS+PHYS_GRID=BLOCK_BUF_NRECS 5   �  @     READ_RESTART%NLCOLS+PHYS_GRID=NLCOLS 5     @     READ_RESTART%NGCOLS+PHYS_GRID=NGCOLS 7   ]  @     READ_RESTART%NCHUNKS+PHYS_GRID=NCHUNKS 6   �  @     READ_RESTART%ENDCHUNK+PPGRID=ENDCHUNK 6   �  @     READ_RESTART%BEGCHUNK+PPGRID=BEGCHUNK 4     @     READ_RESTART%NUMLATS+PMGRID=NUMLATS 0   ]  @     READ_RESTART%NPES+SPMD_DYN=NPES %   �  �       SET_RESTART_FILEPATH .   3  A      SET_RESTART_FILEPATH%LEN_TRIM *   t  =      SET_RESTART_FILEPATH%TRIM ,   �  L   a   SET_RESTART_FILEPATH%RGPATH    �    �   RESTART!COMCTL      H      ITSST    U  H      NSREST    �  H      IRADSW    �  H      IRADLW    -   H      IRADAE    u   H      NREFRQ    �   H      ANNCYC    !  H      NLEND    M!  H      NLRES    �!  H      NLHST    �!  H      LBRNCH    %"  H      AERES    m"  H      OZNCYC    �"  H      SSTCYC    �"  H      ICECYC    E#  H      ADIABATIC    �#  H      FLXAVE    �#  H      IDEAL_PHYS    $  H      NSPLIT    e$  H      IORD    �$  H      JORD    �$  H      KORD    =%  H      USE_ETA    �%  H      AQUA_PLANET    �%  H      DORAMP_SO4    &  H      DORAMP_SCON    ]&  H      FULLGRID     �&  H      PRINT_STEP_COST    �&  H      DOABSEMS    5'  H      DOSW    }'  H      DOLW    �'  H      INDIRECT     (  H      SOM_CONSCHK_FRQ     U(  H      ICE_CONSCHK_FRQ "   �(  �   �   RESTART!COMCTL_R8    )  H      DIVDAMPN    g)  H      PRECC_THRESH    �)  H      PRECL_THRESH    �)  �   �   RESTART!COMLUN    �*  H      NSDS    �*  H      NRG    8+  H      NRG2    �+  H      NCID_INI    �+  H      NCID_OZ    ,  H      NCID_SST    X,  H      NCID_TRC    �,  H      LUHREST    �,  �   �   RESTART!COMHYB    �-  �      HYAI    V.  �      HYAM    �.  �      HYBI    �/  �      HYBM    B0  �      HYBD    �0  �      HYPI    �1  �      HYPM    .2  �      HYPD    �2  H      PS0    3  H      PSR    b3  H      PRSFAC    �3  H      NPRLEV    �3  Y   �   RESTART!COMTFC    K4  H      EPS 