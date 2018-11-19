import maya.cmds as cmds
import maya.mel as mel
import math

# ****************************
# ---------- SETUP ----------
# ****************************

# ----- Create Light -----

# Create ambient light

light = cmds.ambientLight(intensity=1.0)
cmds.ambientLight( light, e=True, ss=True, intensity=0.2, n='lightAmb')
cmds.move( 0, 8, 0 )

# Create directional light
light = cmds.directionalLight(rotation=(45, 30, 15), n='lightDir')
cmds.directionalLight( light, e=True, ss=True, intensity=0.0 )

# Query it
cmds.ambientLight( light, q=True, intensity=True )
cmds.directionalLight( light, q=True, intensity=True )


# ----- Create Transparent Box -----

# Create the groundPlane
cmds.polyCube(w=5, h=2, d=5, sx=1, sy=1, sz=1, ax=(0, 1, 0), name='transparentCube')
cmds.polyNormal() #change normals

# Delete top of cube
cmds.select( 'transparentCube.f[1]' )
cmds.delete()

# Create transparent material for transparent box
cmds.select( 'transparentCube' )
cmds.setAttr( 'lambert1.transparency', 0.9, 0.9, 0.9, type = 'double3' )
cmds.setAttr( 'lambert1.refractions', 1 )
cmds.setAttr( 'lambert1.refractiveIndex', 1.52 )
cmds.setAttr( 'lambert1.color', 0.4, 0.8, 1, type = 'double3' )

# ----- Create Particles -----

count = 0
for i in range( 0, 8 ):
    for j in range( 0, 8 ):
        for k in range( 0, 8 ):
            count=count+1
            result = cmds.polySphere( r=0.15, sx=1, sy=1, name='particle#' )
            cmds.select('particle' + str(count)) 
            cmds.move(-i*0.35+0.75, 3+j*0.35, k*0.35-0.75,'particle' + str(count))


# ****************************
# ------ FOR ANIMATION -------
# ****************************

# ------ FUNCTIONS ------

# Set Keyframes
def setNextKeyParticle( pName, pKeyStart, pTargetAttribute, pValue ):
    
    # keyNext = pKeyStart + pDt
    keyNext = pKeyStart
    
    # clear selection list and select all particles
    cmds.select( clear=True )

    # create animation, set startkeyframe, startvalue=0 at first key frame. Make linear keyframes
    cmds.setKeyframe( pName, time=keyNext, attribute=pTargetAttribute, value=pValue )
    
    # cmds.setKeyframe( pName, time=pEndKey, attribute=pTargetAttribute, value=pValue ) 
    cmds.selectKey( pName, time=( pKeyStart, keyNext ), attribute=pTargetAttribute, keyframe=True )
    cmds.keyTangent( inTangentType='linear', outTangentType='linear' )


# ******************************************************#

# ---------------------- MAIN --------------------------

# ******************************************************#

nrOfParticles = 8*8*8+1
mass = 5

# Playback options
keyFrames = 50
cmds.playbackOptions( playbackSpeed = 0, maxPlaybackSpeed = 1, min = 1, max = 50 )
startTime = cmds.playbackOptions( query = True, minTime = True )
endTime = cmds.playbackOptions( query = True, maxTime = True )
time = startTime
dt = 1 # 24 f/s

# Set first Keyframe for all particles
for i in range ( 1, nrOfParticles ):
    pos = [ cmds.getAttr( 'particle'+str(i)+'.translateX' ),
            cmds.getAttr( 'particle'+str(i)+'.translateY' ),
            cmds.getAttr( 'particle'+str(i)+'.translateZ' ) ]
    
    setNextKeyParticle( 'particle'+str(i), time, 'translateX', pos[0] )
    setNextKeyParticle( 'particle'+str(i), time, 'translateY', pos[1] )
    setNextKeyParticle( 'particle'+str(i), time, 'translateZ', pos[2] )
    
counter = 0
for j in range ( 1, keyFrames ):
    print 'frame: ' + str(j)
    time += dt
    counter = counter + 0.1
    for i in range ( 1, nrOfParticles ):
        
        particlePos = [ cmds.getAttr( 'particle'+str(i)+'.translateX' ),
                        cmds.getAttr( 'particle'+str(i)+'.translateY' ),
                        cmds.getAttr( 'particle'+str(i)+'.translateZ' ) ]
       
        newParticlePosition = [ particlePos[0] + counter, particlePos[1] + counter, particlePos[2] + counter ]

        # Set keyframes
        cmds.select( 'particle'+str(i) )
        setNextKeyParticle( 'particle'+str(i), time, 'translateX', newParticlePosition[0] )
        setNextKeyParticle( 'particle'+str(i), time, 'translateY', newParticlePosition[1] )
        setNextKeyParticle( 'particle'+str(i), time, 'translateZ', newParticlePosition[2] )
