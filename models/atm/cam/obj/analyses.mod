	  �  K   k820309              15.0        ܽW                                                                                                           
       /home2/n02/n02/pappas/src/cam3_sp_intel/models/atm/cam/src/control/analyses.F90 ANALYSES       	       ANALYSES_INI ANALYSES_INT ANALYSES_NUDGE T_A U_A V_A Q_A PS_A S_A                                                    
                                                           
       PLON PLEV PLAT MASTERPROC          @       �                                  
       PCOLS PVER BEGCHUNK ENDCHUNK                      @                              
       SCATTER_FIELD_TO_CHUNK GET_NCOLS_P                      @                              
       GETFIL                      @                              
       NCID_ANALYSES LESS_SURFACE_NUDGING NUDGE_DSE_NOT_T                  � @                              
       R8 SHR_KIND_R8               � @                              '                     #DAY 	   #TOD 
                � D                            	                                � D                             
                          #ESMF_TOD                   � @                              '                    #TYPE    #SEC    #MSEC                 � D                                                            � D                                                           � D                                                            � @                              '                    #CALENDAR    #YEAR    #MONTH    #DAY    #TOD    #JULIANDAY    #DAYOFYEAR                 � D                                  �                      #ESMF_CALENDAR                   � @                              '�                    #TYPE    #DIM    #DIMRUNNINGSUM    #DIY                 � D                                                            � D                                                           p          p            p                                       � D                                        p                   p          p            p                                       � D                                 �                          � D                                 �                          � D                                 �                          � D                                 �                          � D                                         �              #ESMF_TOD                 � D                                                          � D                                                           � @                               '�                   #NSTEP    #STEPSIZE    #STARTDATE    #STOPDATE     #BASEDATE !   #CURRDATE "   #PREVDATE #                � D                                                            � D                                                        #ESMF_TIME                 � D                                         (              #ESMF_DATE                 � D                                          H             #ESMF_DATE                 � D                             !            h             #ESMF_DATE                 � D                             "            �             #ESMF_DATE                 � D                             #            �             #ESMF_DATE    #         @                                  $                   #GETFIL%LEN_TRIM %   #GETFIL%TRIM &   #GETFIL%PRESENT '   #FULPATH (   #LOCFN )   #IFLAG *                 @                           %     LEN_TRIM               @                           &     TRIM               @                           '     PRESENT           
                                 (                    1                                           )                     1           
                                 *                      @@                               +                      @                                 ,                      @                                 -                                                      .            #         @                                   /                    #ANALYSES_INI%LOG 0   #ANALYSES_INI%ABS 1   #ANALYSES_INI%TRIM 2                 @                           0     LOG               @                           1     ABS               @                           2     TRIM #         @                                   3                   #ANALYSES_INT%EPSILON 4   #ANALYSES_INT%LOG 5   #ANALYSES_INT%ABS 6   #ANALYSES_INT%TRIM 7   #ZTODT 8                 @                           4     EPSILON               @                           5     LOG               @                           6     ABS               @                           7     TRIM           
                                 8     
      #         @                                   9                   #ANALYSES_NUDGE%EXP :   #ANALYSES_NUDGE%MAX ;   #ZTODT <   #TAU =   #NCOL >   #NVER ?   #ANALYSIS @   #FIELD A   #TEND B   #L_TEND C                 @                           :     EXP               @                           ;     MAX           
                                 <     
                
                                 =     
                
                                  >                     
                                  ?                    
                                 @                    
 "   p 	         p          5 � p        r ?       p          5 � p        r ?                              
                                 A                    
 #   p 	         p          5 � p        r ?       p          5 � p        r ?                              D                                B                    
 $    p 	         p          5 � p        r ?       p          5 � p        r ?                               D                                 C                     @                                D                   
                &                   &                   &                                                    @                                E                   
                &                   &                   &                                                    @                                F                   
                &                   &                   &                                                    @                                G                   
                &                   &                   &                                                    @                                H                   
                &                   &                                                    @                                I                   
                &                   &                   &                                              �   a      fn#fn      R   b   uapp(ANALYSES    S  @   J  MPISHORTHAND    �  Z   J  PMGRID    �  ]   J  PPGRID    J  c   J  PHYS_GRID    �  G   J  IOFILEMOD    �  s   J  RUNTIME_OPTS    g  O   J  SHR_KIND_MOD '   �  b       ESMF_TIME+ESMF_TIMEMOD /     H   %   ESMF_TIME%DAY+ESMF_TIMEMOD=DAY /   `  ^   %   ESMF_TIME%TOD+ESMF_TIMEMOD=TOD %   �  m      ESMF_TOD+ESMF_TODMOD /   +  H   %   ESMF_TOD%TYPE+ESMF_TODMOD=TYPE -   s  H   %   ESMF_TOD%SEC+ESMF_TODMOD=SEC /   �  H   %   ESMF_TOD%MSEC+ESMF_TODMOD=MSEC '     �       ESMF_DATE+ESMF_DATEMOD 9   �  c   %   ESMF_DATE%CALENDAR+ESMF_DATEMOD=CALENDAR /   	        ESMF_CALENDAR+ESMF_CALENDARMOD 9   �  H   %   ESMF_CALENDAR%TYPE+ESMF_CALENDARMOD=TYPE 7   �  �   %   ESMF_CALENDAR%DIM+ESMF_CALENDARMOD=DIM K   l  �   %   ESMF_CALENDAR%DIMRUNNINGSUM+ESMF_CALENDARMOD=DIMRUNNINGSUM 7   	  H   %   ESMF_CALENDAR%DIY+ESMF_CALENDARMOD=DIY 1   P	  H   %   ESMF_DATE%YEAR+ESMF_DATEMOD=YEAR 3   �	  H   %   ESMF_DATE%MONTH+ESMF_DATEMOD=MONTH /   �	  H   %   ESMF_DATE%DAY+ESMF_DATEMOD=DAY /   (
  ^   %   ESMF_DATE%TOD+ESMF_DATEMOD=TOD ;   �
  H   %   ESMF_DATE%JULIANDAY+ESMF_DATEMOD=JULIANDAY ;   �
  H   %   ESMF_DATE%DAYOFYEAR+ESMF_DATEMOD=DAYOFYEAR -     �       ESMF_TIMEMGR+ESMF_TIMEMGRMOD 9   �  H   %   ESMF_TIMEMGR%NSTEP+ESMF_TIMEMGRMOD=NSTEP ?     _   %   ESMF_TIMEMGR%STEPSIZE+ESMF_TIMEMGRMOD=STEPSIZE A   m  _   %   ESMF_TIMEMGR%STARTDATE+ESMF_TIMEMGRMOD=STARTDATE ?   �  _   %   ESMF_TIMEMGR%STOPDATE+ESMF_TIMEMGRMOD=STOPDATE ?   +  _   %   ESMF_TIMEMGR%BASEDATE+ESMF_TIMEMGRMOD=BASEDATE ?   �  _   %   ESMF_TIMEMGR%CURRDATE+ESMF_TIMEMGRMOD=CURRDATE ?   �  _   %   ESMF_TIMEMGR%PREVDATE+ESMF_TIMEMGRMOD=PREVDATE !   H  �       GETFIL+IOFILEMOD 3   �  A      GETFIL%LEN_TRIM+IOFILEMOD=LEN_TRIM +   .  =      GETFIL%TRIM+IOFILEMOD=TRIM 1   k  @      GETFIL%PRESENT+IOFILEMOD=PRESENT )   �  L   a   GETFIL%FULPATH+IOFILEMOD '   �  L   a   GETFIL%LOCFN+IOFILEMOD '   C  @   a   GETFIL%IFLAG+IOFILEMOD +   �  @       NCID_ANALYSES+RUNTIME_OPTS 2   �  @       LESS_SURFACE_NUDGING+RUNTIME_OPTS -     @       NUDGE_DSE_NOT_T+RUNTIME_OPTS    C  @       NPES+SPMD_DYN    �  �       ANALYSES_INI !     <      ANALYSES_INI%LOG !   J  <      ANALYSES_INI%ABS "   �  =      ANALYSES_INI%TRIM    �  �       ANALYSES_INT %   s  @      ANALYSES_INT%EPSILON !   �  <      ANALYSES_INT%LOG !   �  <      ANALYSES_INT%ABS "   +  =      ANALYSES_INT%TRIM #   h  @   a   ANALYSES_INT%ZTODT    �  �       ANALYSES_NUDGE #   w  <      ANALYSES_NUDGE%EXP #   �  <      ANALYSES_NUDGE%MAX %   �  @   a   ANALYSES_NUDGE%ZTODT #   /  @   a   ANALYSES_NUDGE%TAU $   o  @   a   ANALYSES_NUDGE%NCOL $   �  @   a   ANALYSES_NUDGE%NVER (   �  �   a   ANALYSES_NUDGE%ANALYSIS %   �  �   a   ANALYSES_NUDGE%FIELD $   �  �   a   ANALYSES_NUDGE%TEND &   k  @   a   ANALYSES_NUDGE%L_TEND    �  �       T_A    g  �       U_A    #  �       V_A    �  �       Q_A    �  �       PS_A    ?  �       S_A 