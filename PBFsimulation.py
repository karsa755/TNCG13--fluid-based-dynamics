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
height = 8
width = 8
depth = 8
for i in range( 0, width ):
    for j in range( 0, height ):
        for k in range( 0, depth ):
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


# ******************************************************** #
# -------------------- LINEAR ALGEBRA -------------------- #
# ******************************************************** #

def dotProduct(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def crossProduct(v1, v2):
    return [v1[1]*v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]] 

def addToVec(res, v):
    res[0] = res[0] + v[0]
    res[1] = res[1] + v[1]
    res[2] = res[2] + v[2]

def addVec(v1, v2):
    res = []
    res.append(v1[0] + v2[0])
    res.append(v1[1] + v2[1])
    res.append(v1[2] + v2[2])
    return res
  
def subtractVec(v1, v2):
    res = []
    res.append(v1[0] - v2[0])
    res.append(v1[1] - v2[1])
    res.append(v1[2] - v2[2])
    return res

def scalarMult(v, s) :
    v[0] = s * v[0]
    v[1] = s * v[1]
    v[2] = s * v[2]

def returnScalarMult(v, s) :
    res = []
    res.append(s * v[0])
    res.append(s * v[1])
    res.append(s * v[2])
    return res

def normalize(v) :
    length = getLengthOfVec(v)
    return returnScalarMult(v, 1.0 / length)

def getLengthOfVec(v) :
    return math.sqrt(dotProduct(v,v))

def getDistance(v1, v2) :
    return math.sqrt(math.pow(v1[0]-v2[0], 2) + math.pow(v1[1]-v2[1], 2) + math.pow(v1[2]-v2[2], 2))

# ******************************************************** #
# ----------------------- KERNELS ------------------------ #
# ******************************************************** #

# Pressure kernel (gradient calculations)
def calculateSpikyGradient(h, p1, p2) :
    r = subtractVec(p1, p2)
    r_len = getLengthOfVec(r)

    if r_len <= 0.0 or r_len >= h :
        return [0.0, 0.0, 0.0]

    r = normalize(r)

    spikyConstant = -45.0 / (math.pi * math.pow(h, 6))
    gradX = spikyConstant * math.pow(h-r_len, 2) * r[0]
    gradY = spikyConstant * math.pow(h-r_len, 2) * r[1]
    gradZ = spikyConstant * math.pow(h-r_len, 2) * r[2]

    return [gradX, gradY, gradZ]

# Density kernel
def calculatePoly6(h, p1, p2) :
    r = subtractVec(p1, p2)
    r_len = getLengthOfVec(r)

    return calculatePoly6Scalar(h, r_len)

def calculatePoly6Scalar(h, scalar) :
    if scalar <= 0.0 or scalar >= h :
        return 0.0
    
    poly6Constant = 315.0 / (math.pi * 64.0 * math.pow(h, 9))
    grad = poly6Constant * math.pow(math.pow(h, 2) - math.pow(scalar, 2), 3) 
    return grad

# ******************************************************** #
# ----------- POSITION BASED FLUID SIMULATION ------------ #
# ******************************************************** #

# Find Neighbouring particles and store them in an 2d array
def findNeighbours( posX, posY, posZ, nrOfParticles, radius ):
    neighbours = []
    neighbours.append([0]) # dummy

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
            scalarMult(gradient_j, 1.0/rho_0)
            sum_gradient = sum_gradient + dotProduct(gradient_j, gradient_j)
            addToVec(gradient_i, gradient_j)

            densityResult = densityResult + calculatePoly6(h, n_j, n_i)
        
        sum_gradient = sum_gradient + dotProduct(gradient_i, gradient_i)

        density_constraint = (densityResult / rho_0) - 1.0
        result = (-1.0 * density_constraint) / (sum_gradient + epsilon)
        lambdaArray.append(result)
        
    return lambdaArray

def calculateDeltaPos( posX, posY, posZ, neighbours, nrOfParticles, rho_0, epsilon, lambdas, k, h, n, dQ):
    deltaPosX = []
    deltaPosY = []
    deltaPosZ = []
    deltaPosX.append(0.0) # dummy
    deltaPosY.append(0.0) # dummy
    deltaPosZ.append(0.0) # dummy
    
    for i in range (1, nrOfParticles) :
        lambda_i = lambdas[i]
        pos_i = [posX[i], posY[i], posZ[i]]
        resultX = 0.0
        resultY = 0.0
        resultZ = 0.0
        for j in range (1, len(neighbours[i])) :
            lambda_j = lambdas[neighbours[i][j]]
            lambda_ij = lambda_j + lambda_i

            pos_j = [posX[neighbours[i][j]], posY[neighbours[i][j]], posZ[neighbours[i][j]]]
            corrScore = computeCorrectionScore(k, h, n, dQ, pos_j, pos_i)
            spiky = calculateSpikyGradient(h, pos_j, pos_i)

            resultX = resultX + (lambda_ij + corrScore) * spiky[0] 
            resultY = resultY + (lambda_ij + corrScore) * spiky[1] 
            resultZ = resultZ + (lambda_ij + corrScore) * spiky[2] 
        
        resultX = resultX * 1.0/rho_0
        resultY = resultY * 1.0/rho_0
        resultZ = resultZ * 1.0/rho_0
        deltaPosX.append(resultX)
        deltaPosY.append(resultY)
        deltaPosZ.append(resultZ)

    return [deltaPosX, deltaPosY, deltaPosZ] 

def computeCorrectionScore(k, h, n, dQ, p1, p2) :
    constraint = calculatePoly6(h, p1, p2) / calculatePoly6Scalar(h, dQ * h)
    return -k * math.pow(constraint, n)

def calculateParticleCollisionResponse(posX, posY, posZ, vX, vY, vZ, particleRadius, neighbours) :
    offset = 0.01

    for i in range (1, nrOfParticles) :
        pos_i = [posX[i], posY[i], posZ[i]]
        for j in range (1, len(neighbours[i])) :
            pos_j = [posX[neighbours[i][j]], posY[neighbours[i][j]], posZ[neighbours[i][j]]]
            distance = getDistance(pos_i, pos_j)
            if distance <= 2.0 * particleRadius and distance > 0.0 : # COLLISION!
                v_i = [vX[i], vY[i], vZ[i]]
                v_j = [vX[neighbours[i][j]], vY[neighbours[i][j]], vZ[neighbours[i][j]]]

                vecBetweenParticles = [pos_j[0] - pos_i[0], pos_j[1] - pos_i[1], pos_j[2] - pos_i[2]]
                collisionPoint = addVec(pos_i, vecBetweenParticles)

                vecBetweenParticles = normalize(vecBetweenParticles)

                negVec = returnScalarMult(returnScalarMult(vecBetweenParticles, -1.0), offset)
                posX[i] = posX[i] + negVec[0]
                posY[i] = posY[i] + negVec[1]
                posZ[i] = posZ[i] + negVec[2]
                
                posX[neighbours[i][j]] = posX[neighbours[i][j]] + offset * vecBetweenParticles[0]
                posY[neighbours[i][j]] = posY[neighbours[i][j]] + offset * vecBetweenParticles[1]
                posZ[neighbours[i][j]] = posZ[neighbours[i][j]] + offset * vecBetweenParticles[2]

                newVels = calculateNewVelocities(pos_i, pos_j, v_i, v_j)
                vX[i] = newVels[0][0]
                vY[i] = newVels[0][1]
                vZ[i] = newVels[0][2]
                
                vX[neighbours[i][j]] = newVels[1][0]
                vY[neighbours[i][j]] = newVels[1][1]
                vZ[neighbours[i][j]] = newVels[1][2]

    return [posX, posY, posZ, vX, vY, vZ]

def calculateNewVelocities(pos_i, pos_j, v_i, v_j) :
    pos_ij = subtractVec(pos_i, pos_j)
    pos_ji = subtractVec(pos_j, pos_i)

    nv1 = projectV1V2(v_j, pos_ji)
    nv1 = subtractVec(nv1, projectV1V2(v_i, pos_ij))

    nv2 = projectV1V2(v_i, pos_ji)
    nv2 = subtractVec(nv2, projectV1V2(v_j, pos_ij))

    scaleVel = 0.25
    scalarMult(nv1, scaleVel)
    scalarMult(nv2, scaleVel)

    return [nv1, nv2]

def projectV1V2(v1, v2) :
    c_1 = dotProduct(v1, v2)
    c_2 = dotProduct(v2, v2)

    finalV = returnScalarMult(v1, c_1)
    scalarMult(finalV, 1.0/c_2)

    return finalV    



def calculateCollisionResponse(posX, posY, posZ, vX, vY, vZ, particleRadius) :
    # Bounding condition for the transparent box
    xMin = -2.5 + particleRadius
    xMax = 2.5 - particleRadius
    yMin = -0.5
    zMin = -2.5 + particleRadius
    zMax = 2.5 - particleRadius

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
    vorticity.append([0.0, 0.0, 0.0]) # dummy

    for i in range (1, nrOfParticles) :
        result = [0.0, 0.0, 0.0]
        v_i = [vX[i], vY[i], vZ[i]]
        pos_i = [posX[i], posY[i], posZ[i]]
        for j in range (1, len(neighbours[i])) :
            v_j = [vX[neighbours[i][j]], vY[neighbours[i][j]], vZ[neighbours[i][j]]]
            pos_j = [posX[neighbours[i][j]], posY[neighbours[i][j]], posZ[neighbours[i][j]]]
            v_ij = [v_j[0] - v_i[0], v_j[1] - v_i[1], v_j[2] - v_i[2]]
            spiky = calculateSpikyGradient(h, pos_i, pos_j)
            addToVec(result, crossProduct(v_ij, spiky))

        vorticity.append(result)
        
    return vorticity

def computeCorrectiveForce(h, posX, posY, posZ, vorticity, vortEps) :
    corrForce = []
    corrForce.append([0.0, 0.0, 0.0]) # dummy

    for i in range (1, nrOfParticles) :
        result = [0.0, 0.0, 0.0]

        pos_i = [posX[i], posY[i], posZ[i]]
        for j in range (1, len(neighbours[i])) :
            pos_j = [posX[neighbours[i][j]], posY[neighbours[i][j]], posZ[neighbours[i][j]]]
            spiky = calculateSpikyGradient(h, pos_i, pos_j)
            vortLength = getLengthOfVec(vorticity[neighbours[i][j]])
            scalarMult(spiky, vortLength)
            addToVec(result, spiky)
        
        resultLength = getLengthOfVec(result)
        
        if resultLength <= 0.0001 :
            corrForce.append([0.0, 0.0, 0.0])
            continue

        result = normalize(result)

        cross = crossProduct(result, vorticity[i])
        scalarMult(cross, vortEps)
        corrForce.append(cross)

    return corrForce
            

def computeXSPHViscosity(c, h, posX, posY, posZ, vX, vY, vZ, neighbours) : 
    viscosity = []
    viscosity.append([0.0, 0.0, 0.0]) # dummy

    for i in range (1, nrOfParticles) :
        result = [0.0, 0.0, 0.0]

        v_i = [vX[i], vY[i], vZ[i]]
        pos_i = [posX[i], posY[i], posZ[i]]
        for j in range (1, len(neighbours[i])) :
            v_j = [vX[neighbours[i][j]], vY[neighbours[i][j]], vZ[neighbours[i][j]]]
            pos_j = [posX[neighbours[i][j]], posY[neighbours[i][j]], posZ[neighbours[i][j]]]
            polyKernel = calculatePoly6(h, pos_j, pos_i)
            v_ij = [v_j[0] - v_i[0], v_j[1] - v_i[1], v_j[2] - v_i[2]]
            scalarMult(v_ij, polyKernel)

            addToVec(result, v_ij)
        
        scalarMult(result, c)
        viscosity.append(result)
    
    return viscosity



# ****************************************************** #
# ---------------------- MAIN -------------------------- #
# ****************************************************** #

maxIterations = 10

nrOfParticles = height * width * depth + 1
mass = 1.0
# Velocities
vX = [0] * nrOfParticles
vY = [0] * nrOfParticles
vZ = [0] * nrOfParticles
# Predicted positions
ppX = [0] * nrOfParticles
ppY = [0] * nrOfParticles
ppZ = [0] * nrOfParticles

gravity = -9.82
dt = 0.016
rho_0 = 1000.0
epsilon = 200.0
correctionK = 0.001
correctionN = 4.0
correctionDeltaQ = 0.3
vorticityEps = 0.25
XSPHC = 0.001
h = 0.4

# Playback options
keyFrames = 150
cmds.playbackOptions( playbackSpeed = 0, maxPlaybackSpeed = 1, min = 1, max = 150 )
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

        # Impose boundary constraint
        collision = calculateCollisionResponse(ppX[i], ppY[i], ppZ[i], vX[i], vY[i], vZ[i], particleRadius)
        ppX[i] = collision[0]
        ppY[i] = collision[1]
        ppZ[i] = collision[2]
        vX[i] = collision[3]
        vY[i] = collision[4]
        vZ[i] = collision[5]

    neighbours = findNeighbours(ppX, ppY, ppZ, nrOfParticles, h)
    
    iterations = 0
    while iterations < maxIterations :
        lambdas = calculateLambda( ppX, ppY, ppZ, neighbours, nrOfParticles, rho_0, epsilon, h)     
        deltaPos = calculateDeltaPos( ppX, ppY, ppZ, neighbours, nrOfParticles, rho_0, epsilon, lambdas, correctionK, h, correctionN, correctionDeltaQ)

        for i in range (1, nrOfParticles) :
            ppX[i] = ppX[i] + deltaPos[0][i]
            ppY[i] = ppY[i] + deltaPos[1][i]
            ppZ[i] = ppZ[i] + deltaPos[2][i]

        iterations = iterations + 1

    particleCollision = calculateParticleCollisionResponse(ppX, ppY, ppZ, vX, vY, vZ, particleRadius, neighbours)
    ppX = particleCollision[0]
    ppY = particleCollision[1]
    ppZ = particleCollision[2]
    vX = particleCollision[3]
    vY = particleCollision[4]
    vZ = particleCollision[5]

    for i in range (1, nrOfParticles) :           
        collision = calculateCollisionResponse(ppX[i], ppY[i], ppZ[i], vX[i], vY[i], vZ[i], particleRadius)
        ppX[i] = collision[0]
        ppY[i] = collision[1]
        ppZ[i] = collision[2]
        vX[i] = collision[3]
        vY[i] = collision[4]
        vZ[i] = collision[5]

    for i in range (1, nrOfParticles) :
        pos = [ cmds.getAttr( 'particle'+str(i)+'.translateX' ),
                cmds.getAttr( 'particle'+str(i)+'.translateY' ),
                cmds.getAttr( 'particle'+str(i)+'.translateZ' ) ]

        vX[i] = (1.0/dt) * (ppX[i] - pos[0])
        vY[i] = (1.0/dt) * (ppY[i] - pos[1])
        vZ[i] = (1.0/dt) * (ppZ[i] - pos[2])
    
    vorticity = calculateVorticityConfinement(h, ppX, ppY, ppZ, vX, vY, vZ, neighbours)
    corrForce = computeCorrectiveForce(h, ppX, ppY, ppZ, vorticity, vorticityEps)
    viscosity = computeXSPHViscosity(XSPHC, h, ppX, ppY, ppZ, vX, vY, vZ, neighbours)

    for i in range (1, nrOfParticles) :
        vX[i] = vX[i] + viscosity[i][0] + dt*corrForce[i][0] 
        vY[i] = vY[i] + viscosity[i][1] + dt*corrForce[i][1] 
        vZ[i] = vZ[i] + viscosity[i][2] + dt*corrForce[i][2] 
    

    for i in range (1, nrOfParticles) :
        cmds.select( 'particle'+str(i) )
        setNextKeyParticle( 'particle'+str(i), time, 'translateX', ppX[i] )
        setNextKeyParticle( 'particle'+str(i), time, 'translateY', ppY[i] )
        setNextKeyParticle( 'particle'+str(i), time, 'translateZ', ppZ[i] )
