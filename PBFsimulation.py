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

# Find Neighbouring particles and store them in an 2d array
def findNeighbours( posX, posY, posZ, nrOfParticles ):
    radius = 0.5
    neighbours = []
    neighbours.append([0])

    for i in range ( 1, nrOfParticles ): 
        closestNeighbours = []

        for j in range (1, nrOfParticles ):
            if i == j :
                continue
            arrayDistance = [ posX[i] - posX[j], posY[i] - posY[j], posZ[i] - posZ[j] ]
            distance = math.sqrt(math.pow(arrayDistance[0], 2) + math.pow(arrayDistance[1], 2) + math.pow(arrayDistance[2], 2))

            if(distance < radius) :
                closestNeighbours.append(j)
        
        neighbours.append(closestNeighbours)
        
    return neighbours

# Calculate lamda (particle constraint) to enforce incompressibility
def calculateLambda( posX, posY, posZ, neighbours, nrOfParticles) :
    p_constraint = 0
    
    return p_constraint


# Pressure kernel (gradient calculations)
def calculateSpikyGradient(h, p1, p2) :
    deltaPosX = math.abs(p2[0] - p1[0])
    deltaPosY = math.abs(p2[1] - p1[1])
    deltaPosZ = math.abs(p2[2] - p1[2])

    spikyConstant = 15 / (math.pi * math.pow(h, 6))
    gradX = spikyConstant * math.pow(h-deltaPosX, 3) if deltaPosX >= 0 & deltaPosX <= h else 0
    gradY = spikyConstant * math.pow(h-deltaPosY, 3) if deltaPosY >= 0 & deltaPosY <= h else 0
    gradX = spikyConstant * math.pow(h-deltaPosZ, 3) if deltaPosZ >= 0 & deltaPosZ <= h else 0

    return [gradX, gradY, gradZ]

# Density kernel
def calculatePoly6(h, p1, p2) :
    deltaPosX = math.abs(p2[0] - p1[0])
    deltaPosY = math.abs(p2[1] - p1[1])
    deltaPosZ = math.abs(p2[2] - p1[2])

    poly6Constant = 315 / (math.pi * 64 * math.pow(h, 9))
    gradX = poly6Constant * math.pow(math.pow(h, 2) - math.pow(deltaPosX, 2), 3) if deltaPosX >= 0 & deltaPosX <= h else 0
    gradY = poly6Constant * math.pow(math.pow(h, 2) - math.pow(deltaPosY, 2), 3) if deltaPosY >= 0 & deltaPosY <= h else 0
    gradZ = poly6Constant * math.pow(math.pow(h, 2) - math.pow(deltaPosZ, 2), 3) if deltaPosZ >= 0 & deltaPosZ <= h else 0

    return [gradX, gradY, gradZ]

# ******************************************************#

# ---------------------- MAIN --------------------------

# ******************************************************#

maxIterations = 10

nrOfParticles = 8*8*8+1
mass = 5
# Velocities
vX = [0] * nrOfParticles
vY = [0] * nrOfParticles
vZ = [0] * nrOfParticles
# Predicted positions
ppX = [0] * nrOfParticles
ppY = [0] * nrOfParticles
ppZ = [0] * nrOfParticles
gravity = -9.82
dt = 0.1


# Playback options
keyFrames = 5
cmds.playbackOptions( playbackSpeed = 0, maxPlaybackSpeed = 1, min = 1, max = 5 )
startTime = cmds.playbackOptions( query = True, minTime = True )
endTime = cmds.playbackOptions( query = True, maxTime = True )
time = startTime
keyStep = 1


# Set first Keyframe for all particles
for i in range ( 1, nrOfParticles ):
    pos = [ cmds.getAttr( 'particle'+str(i)+'.translateX' ),
            cmds.getAttr( 'particle'+str(i)+'.translateY' ),
            cmds.getAttr( 'particle'+str(i)+'.translateZ' ) ]

    ppX[i] = pos[0]
    ppY[i] = pos[1]
    ppZ[i] = pos[2]
    
    setNextKeyParticle( 'particle'+str(i), time, 'translateX', pos[0] )
    setNextKeyParticle( 'particle'+str(i), time, 'translateY', pos[1] )
    setNextKeyParticle( 'particle'+str(i), time, 'translateZ', pos[2] )
     

for j in range ( 1, keyFrames ):
    print 'frame: ' + str(j)
    time += keyStep
    
    for i in range (1, nrOfParticles) :
        # Apply external forces to velocity
        vY[i] = vY[i] + dt*gravity
        # Apply velocity to position estimate
        ppY[i] = ppY[i] + dt * vY[i]

    neighbours = findNeighbours(ppX, ppY, ppZ, nrOfParticles)
    iterations = 0
    while iterations < maxIterations :
        
        
        
        iterations = iterations + 1

    for i in range (1, nrOfParticles) :
        cmds.select( 'particle'+str(i) )
        setNextKeyParticle( 'particle'+str(i), time, 'translateX', ppX[i] )
        setNextKeyParticle( 'particle'+str(i), time, 'translateY', ppY[i] )
        setNextKeyParticle( 'particle'+str(i), time, 'translateZ', ppZ[i] )




 
