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
particleRadius = 0.15
count = 0
for i in range( 0, 8 ):
    for j in range( 0, 8 ):
        for k in range( 0, 8 ):
            count=count+1
            result = cmds.polySphere( r=particleRadius, sx=1, sy=1, name='particle#' )
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
def calculateLambda( posX, posY, posZ, neighbours, nrOfParticles, rho_0, epsilon, h) :

    lambdaArray = []
    lambdaArray.append(0) #dummy

    for i in range (1, nrOfParticles) :
        gradient_i = [0, 0, 0]
        sum_gradient = 0 
        n_i = [posX[i], posY[i], posZ[i]]

        densityResult = 0

        for j in range (1, len(neighbours[i])) :
            n_j = [posX[neighbours[i][j]], posY[neighbours[i][j]], posZ[neighbours[i][j]]]
            gradient_j = calculateSpikyGradient(h, n_j, n_i)
            gradient_j = [gradient_j[0] * rho_0, gradient_j[1] * rho_0, gradient_j[2] * rho_0] 
            sum_gradient = sum_gradient + dotProduct(gradient_j, gradient_j)
            addToVec(gradient_i, gradient_j)

            densityResult = densityResult + calculatePoly6(h, n_j, n_i)
        
        sum_gradient = sum_gradient + dotProduct(gradient_i, gradient_i)

        density_constraint = (densityResult / rho_0) - 1
        result = -density_constraint / (sum_gradient + epsilon)
        lambdaArray.append(result)
        
    return lambdaArray


def dotProduct(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def crossProduct(v1, v2):
    return [v1[1]*v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]] 

def addToVec(res, v):
    res[0] = res[0] + v[0]
    res[1] = res[1] + v[1]
    res[2] = res[2] + v[2]
  
    

# Pressure kernel (gradient calculations)
def calculateSpikyGradient(h, p1, p2) :
    deltaPosX = p2[0] - p1[0]
    deltaPosY = p2[1] - p1[1]
    deltaPosZ = p2[2] - p1[2]

    spikyConstant = 15 / (math.pi * math.pow(h, 6))
    gradX = spikyConstant * math.pow(h-deltaPosX, 3) if deltaPosX >= 0.0 and deltaPosX <= h else 0.0
    gradY = spikyConstant * math.pow(h-deltaPosY, 3) if deltaPosY >= 0.0 and deltaPosY <= h else 0.0
    gradZ = spikyConstant * math.pow(h-deltaPosZ, 3) if deltaPosZ >= 0.0 and deltaPosZ <= h else 0.0

    return [gradX, gradY, gradZ]

# Density kernel
def calculatePoly6(h, p1, p2) :
    deltaPosX = p2[0] - p1[0]
    deltaPosY = p2[1] - p1[1]
    deltaPosZ = p2[2] - p1[2]

    length = math.sqrt(math.pow(deltaPosX, 2) + math.pow(deltaPosY, 2) + math.pow(deltaPosZ, 2))

    return calculatePoly6Scalar(h, length)

def calculatePoly6Scalar(h, scalar) :
    poly6Constant = 315 / (math.pi * 64 * math.pow(h, 9))
    grad = poly6Constant * math.pow(math.pow(h, 2) - math.pow(scalar, 2), 3) if scalar >= 0 and scalar <= h else 0
    return grad

def calculateDeltaPos( posX, posY, posZ, neighbours, nrOfParticles, rho_0, epsilon, lambdas, k, h, n, dQ):
    deltaPosX = []
    deltaPosY = []
    deltaPosZ = []

    deltaPosX.append(0.0)
    deltaPosY.append(0.0)
    deltaPosZ.append(0.0)
    
    for i in range (1, nrOfParticles) :
        lambda_i = lambdas[i]
        resultX = 0.0
        resultY = 0.0
        resultZ = 0.0
        for j in range (1, len(neighbours[i])) :
            lambda_j = lambdas[neighbours[i][j]]
            lambda_ij = lambda_j + lambda_i

            pos_i = [posX[i], posY[i], posZ[i]]
            pos_j = [posX[neighbours[i][j]], posY[neighbours[i][j]], posZ[neighbours[i][j]]]
            corrScore = computeCorrectionScore(k, h, n, dQ, pos_i, pos_j)
            spiky = calculateSpikyGradient(h, pos_j, pos_i)

            resultX = (lambda_ij + corrScore) * spiky[0] 
            resultY = (lambda_ij + corrScore) * spiky[1] 
            resultZ = (lambda_ij + corrScore) * spiky[2] 
        
        resultX = resultX * rho_0
        resultY = resultY * rho_0
        resultZ = resultZ * rho_0
        deltaPosX.append(resultX)
        deltaPosY.append(resultY)
        deltaPosZ.append(resultZ)

    return [deltaPosX, deltaPosY, deltaPosZ] 

def computeCorrectionScore(k, h, n, dQ, p1, p2) :
    constraint = calculatePoly6(h, p1, p2) / calculatePoly6Scalar(h, dQ * h)
    return -k * math.pow(constraint, n)

# MAYBE TODO: Add collisions between particles as well 
def calculateCollisionResponse(posX, posY, posZ, dX, dY, dZ, vX, vY, vZ, particleRadius) :
    # Bounding condition for the transparent box
    xMin = -2.5 + particleRadius
    xMax = 2.5 - particleRadius
    yMin = -0.5
    zMin = -2.5 + particleRadius
    zMax = 2.5 - particleRadius

    posX = posX + dX
    posY = posY + dY
    posZ = posZ + dZ

    if posX < xMin :
        posX = xMin
        vX = 0.0
    elif posX > xMax :
        posX = xMax
        vX = 0.0

    if posY < yMin :
        posY = yMin
        vY = 0.0

    if posZ < zMin :
        posZ = zMin
        vZ = 0.0
    elif posZ > zMax :
        posZ = zMax
        vZ = 0.0

    return [posX, posY, posZ, vX, vY, vZ]
    
def calculateVorticityConfinement(h, posX, posY, posZ, vX, vY, vZ, neighbours) :
    vorticity = []
    vorticity.append([0.0, 0.0, 0.0])

    for i in range (1, nrOfParticles) :
        result = [0.0, 0.0, 0.0]
        v_i = [vX[i], vY[i], vZ[i]]
        pos_i = [posX[i], posY[i], posZ[i]]
        for j in range (1, len(neighbours[i])) :
            v_j = [vX[neighbours[i][j]], vY[neighbours[i][j]], vZ[neighbours[i][j]]]
            pos_j = [posX[neighbours[i][j]], posY[neighbours[i][j]], posZ[neighbours[i][j]]]
            v_ij = [v_j[0] - v_i[0], v_j[1] - v_i[1], v_j[2] - v_i[2]]
            spiky = calculateSpikyGradient(h, pos_j, pos_i)
            addToVec(result, crossProduct(v_ij, spiky))

        vorticity.append(result)
        
    return vorticity


# ******************************************************#

# ---------------------- MAIN --------------------------

# ******************************************************#

maxIterations = 2

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
rho_0 = 1.0
epsilon = 200
correctionK = 0.001
correctionN = 4.0
correctionDeltaQ = 0.3
vorticityEps = 1.0
XSPHC = 0.001
h = 1.2

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
        lambdas = calculateLambda( ppX, ppY, ppZ, neighbours, nrOfParticles, rho_0, epsilon, h)     
        deltaPos = calculateDeltaPos( ppX, ppY, ppZ, neighbours, nrOfParticles, rho_0, epsilon, lambdas, correctionK, h, correctionN, correctionDeltaQ)

        for i in range (1, nrOfParticles) :
            collision = calculateCollisionResponse(ppX[i], ppY[i], ppZ[i], deltaPos[0][i], deltaPos[1][i], deltaPos[2][i], vX[i], vY[i], vZ[i], particleRadius)
            ppX[i] = collision[0]
            ppY[i] = collision[1]
            ppZ[i] = collision[2]
            vX[i] = collision[3]
            vY[i] = collision[4]
            vZ[i] = collision[5]

        iterations = iterations + 1

    for i in range (1, nrOfParticles) :
        pos = [ cmds.getAttr( 'particle'+str(i)+'.translateX' ),
                cmds.getAttr( 'particle'+str(i)+'.translateY' ),
                cmds.getAttr( 'particle'+str(i)+'.translateZ' ) ]

        vX[i] = (1.0/dt) * (ppX[i] - pos[0])
        vY[i] = (1.0/dt) * (ppY[i] - pos[1])
        vZ[i] = (1.0/dt) * (ppZ[i] - pos[2])
    
    vorticity = calculateVorticityConfinement(h, ppX, ppY, ppZ, vX, vY, vZ, neighbours)

    for i in range (1, nrOfParticles) :
        cmds.select( 'particle'+str(i) )
        setNextKeyParticle( 'particle'+str(i), time, 'translateX', ppX[i] )
        setNextKeyParticle( 'particle'+str(i), time, 'translateY', ppY[i] )
        setNextKeyParticle( 'particle'+str(i), time, 'translateZ', ppZ[i] )




 
