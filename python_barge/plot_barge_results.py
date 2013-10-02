"""
This program is distributed AS IS, WITHOUT WARRANTY, and NO PROMISE FOR SUPPORT.
USE AT YOUR OWN RISK
"""

if __name__ == '__main__':  
   import matplotlib.pyplot as plt
   import math
   import numpy as np
   import MAP
   import MAPSupport as PyLibMAP

   map_fair1, map_anch1, map_fair2, map_anch2 = ([] for i in range(4))   
   time, fair1, anch1, fair2, anch2, surge, sway, heave, roll, pitch, yaw = ([] for i in range(11))
   
   R        = np.matrix([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) # converts body frame to inertial (global) frame
   r        = np.matrix('0.0;0.0;0.0')                               # platform displacement 
   ri       = np.matrix('0.0;0.0;0.0')                               # fairlead position in body frame
   position = np.matrix('0.0;0.0;0.0')                               # fairlead position in global frame


   with open('test20.out') as f:
      for l in f:
         line = l.strip().split("\t")
         
         time.append ( float(line[0])  )
         fair1.append( float(line[1])  )
         anch1.append( float(line[2])  )
         fair2.append( float(line[3])  )
         anch2.append( float(line[4])  )
         surge.append( float(line[17]) )
         sway.append ( float(line[18]) )
         heave.append( float(line[19]) )
         roll.append ( float(line[20]) )
         pitch.append( float(line[21]) )
         yaw.append  ( float(line[22]) )

   InitIn = MAP.MAP_InitInputType( )
   InitOut  = MAP.MAP_InitOutputType( )

   PyLibMAP.read_input_file( "NREL_barge.map" , InitIn )
   InitIn.SetDepth     ( -150 );
   InitIn.SetGravity   ( 9.81 );
   InitIn.SetSeaDensity( 1020 );
    
   d   = MAP.MAP_OtherStateType( )
   u   = MAP.MAP_InputType( )
   p   = MAP.MAP_ParameterType( )
   z   = MAP.MAP_ConstraintStateType( )
   y   = MAP.MAP_OutputType( )
   msg = MAP.MAP_Message( )
   err = MAP.MAP_ErrStat( )
   
   MAP.MSQS_Init( InitIn, u, p, None, None, z, d, y, None, InitOut, err, msg )
   if err.error_status( ) != 0 :        
      print msg.status( )        

   uval = u.get( )

   for i in range( len(time) ) :
      r  = [ [surge[i]], [sway[i]], [heave[i]] ] # platform position in global frame
      ri = [ [20.0], [20.0], [-4.0] ]            # mooring line position in body frame
      
      phi = roll[i]*math.pi/180
      the = pitch[i]*math.pi/180
      psi = yaw[i]*math.pi/180
   
      cphi = math.cos( phi );   sphi = math.sin( phi );
      cthe = math.cos( the );   sthe = math.sin( the );
      cpsi = math.cos( psi );   spsi = math.sin( psi );
   
      # first column       # second column                       # third column
      R[0,0] = cpsi*cthe;  R[0,1] = cpsi*sthe*sphi - spsi*cphi;  R[0,2] = cpsi*sthe*cphi + spsi*sphi
      R[1,0] = spsi*cthe;  R[1,1] = spsi*sthe*sphi + cpsi*cphi;  R[1,2] = spsi*sthe*cphi - cpsi*sphi
      R[2,0] = -sthe;      R[2,1] = cthe*sphi;                   R[2,2] = cthe*cphi
        
      position = r + R*ri

      uval[0] = position[0]
      uval[1] = position[1]
      uval[2] = position[2]
      uval[3] = position[0]
      uval[4] = position[1]
      uval[5] = position[2]
      u.set(uval)

      MAP.MSQS_UpdateStates( time[i], 0, u, p, None, None, z, d, err, msg )
      if err.error_status( ) !=0 :
         print msg.status( )
    
      MAP.MSQS_CalcOutput( time[i], u, p, None, None, z, d, y, err, msg )
      if err.error_status( ) != 0 :
         print msg.status( )

   # Destroy objects
   MAP.MSQS_End( u, p, None, None, z, d, y, err, msg )
   if err.error_status( ) != 0 :
      print msg.status( )

   with open( 'map.out' ) as f:
      for _ in xrange( 5 ):
         next( f )
      for line in f:
         words = line.split( )
         try :
            map_fair1.append( float(words[1]) )
            map_anch1.append( float(words[2]) )
            map_fair2.append( float(words[3]) )
            map_anch2.append( float(words[4]) )
         except IndexError :
            print 'Blank spot'
         except ValueError :
            print 'Not a numeric value'

plt.plot( time,fair1,'--',color='b',lw=4,alpha=0.5,label='FAST 7.2 (line 1)' )
plt.plot( time,map_fair1,color='b',label='MAP (line 1)' )
plt.plot( time,fair2,'--',color='r',lw=4,alpha=0.5,label='FAST 7.2 (line 2)' )
plt.plot( time,map_fair2,color='r',label='MAP (line 2)' )
plt.legend( )
plt.xlim([0.0, 60.0])
plt.ylabel('Fairlead Force [kN]')
plt.xlabel('Time [s]')
plt.show( )
