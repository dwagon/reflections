#!/usr/bin/env python
#
# Script to play around with genetic algorithms
#
import os
import sys
import tempfile
import time
import math
from pygene.population import Population
from pygene.organism import MendelOrganism
from pygene.gene import FloatGene
import Image
import ImageChops

worldWidth = 20
worldHeight = 20
worldDepth = 20
imageWidth = 64
imageHeight = 64
gendir = "gendir-%d" % os.getpid()
numreflectors = 10

ri = Image.open("reference.png")
reference_image = ri.resize((imageWidth, imageHeight))

generation = 0


###############################################################################
###############################################################################
###############################################################################
class colourGene(FloatGene):
    mutProb = 0.05
    muyAmt = 0.01
    randMin = 0.0
    randMax = 1.0

    ###########################################################################
    def __repr__(self):
        return "C%0.2f" % self.value


###############################################################################
###############################################################################
###############################################################################
class heightGene(FloatGene):
    mutProb = 0.05
    randMin = 0.0
    randMax = worldHeight

    ###########################################################################
    def __repr__(self):
        return "H%0.2f" % self.value


###############################################################################
###############################################################################
###############################################################################
class depthGene(FloatGene):
    mutProb = 0.05
    randMin = 0.0
    randMax = worldDepth

    ###########################################################################
    def __repr__(self):
        return "D%0.2f" % self.value


###############################################################################
###############################################################################
###############################################################################
class radiusGene(FloatGene):
    mutProb = 0.05
    randMin = 0.5
    randMax = 5.0

    ###########################################################################
    def __repr__(self):
        return "R%0.2f" % self.value


###############################################################################
###############################################################################
###############################################################################
class widthGene(FloatGene):
    mutProb = 0.05
    randMin = 0.0
    randMax = worldWidth

    ###########################################################################
    def __repr__(self):
        return "W%f" % self.value


###############################################################################
###############################################################################
###############################################################################
class reflectBox(MendelOrganism):
    genome = {}
    for i in range(numreflectors):
        genome["x_%d" % i] = widthGene
        genome["y_%d" % i] = heightGene
        genome["z_%d" % i] = depthGene
        genome["r_%d" % i] = radiusGene
        genome["cr_%d" % i] = colourGene
        genome["cg_%d" % i] = colourGene
        genome["cb_%d" % i] = colourGene

    ###########################################################################
    def __repr__(self):
        delta = 0.01
        wwd = worldWidth - delta
        wdd = worldDepth - delta
        str = ""
        str += '#include "colors.inc"\n'
        str += '#include "finish.inc"\n'
        str += 'camera { location <-0.5, -0.5, -0.5> look_at  <9, 9, 9> angle 90 }\n'
        str += 'global_settings { ambient_light White }\n'
        str += 'light_source { <%s, %s, %s> color White }\n' % (-1+delta, -1+delta, -1+delta)
        str += 'light_source { <%s, %s, %s> color White }\n' % (wwd, -1+delta, -1+delta)
        str += 'light_source { <%s, %s, %s> color White }\n' % (wwd, -1+delta, worldDepth+delta)
        str += 'light_source { <%s, %s, %s> color White }\n' % (wwd, worldHeight-delta, wdd)
        str += 'plane { <0, 1, 0>, %s pigment { Black }  }\n' % (worldHeight+1)    # top
        str += 'plane { <0, 1, 0>, -1 pigment { White }  }\n'   # bottom
        str += 'plane { <1, 0, 0>, -1 pigment { Black }  }\n'  # left
        str += 'plane { <1, 0, 0>, %s pigment { White }  }\n' % (worldWidth+1)   # right
        str += 'plane { <0, 0, 1>, -1 pigment { Black }  }\n'  # front
        str += 'plane { <0, 0, 1>, %s pigment { Yellow }  }\n' % (worldDepth+1)   # back
        for o in range(numreflectors):
            colour = "rgb <%f, %f, %f>" % (self["cr_%d" % o], self["cg_%d" % o], self["cb_%d" % o])
            sphere = "sphere { <%f, %f, %f>, %f texture { pigment { color %s } }}\n" % \
                (self["x_%d" % o], self["y_%d" % o], self["z_%d" % o], self["r_%d" % o], colour)
            str += sphere
        return str

    ###########################################################################
    def generateScene(self):
        f = tempfile.NamedTemporaryFile(prefix="%s-" % generation, delete=False, suffix='.pov')
        f.write(str(self))
        f.close()
        return f.name

    ###########################################################################
    def renderScene(self, fname):
        pngf = fname.replace('.pov', '.png')
        f = os.popen('povray -D +H%s +W%s +O%s %s 2>&1' % (imageHeight, imageWidth, pngf, fname))
        log = f.read()
        x = f.close()
        if x:
            sys.stderr.write(log)
        return pngf

    ###########################################################################
    def analyseScene(self, fname):
        # Lower fitness is better
        image = Image.open(fname)
        diffimage = ImageChops.difference(reference_image, image)
        h = diffimage.histogram()
        sq = (value*(idx**2) for idx, value in enumerate(h))
        sum_of_squares = sum(sq)
        diff = math.sqrt(sum_of_squares/float(image.size[0] * image.size[1]))
        diff = diff*diff   # Make differences big
        return diff

    ###########################################################################
    def fitness(self):
        povfile = self.generateScene()
        pngfile = self.renderScene(povfile)
        fitness = self.analyseScene(pngfile)
        self.cleanup(povfile, pngfile)
        for o in range(numreflectors):
            fitness += (100.0*self["r_%d" % o])
        return fitness

    ###########################################################################
    def cleanup(self, *args):
        for f in args:
            os.unlink(f)


###############################################################################
###############################################################################
###############################################################################
class reflectorPopulation(Population):
    initPopulation = 20
    species = reflectBox
    childCull = 2
    childCount = 2
    mutants = 0.20
    numNewOrganisms = 2
    incest = 2


###############################################################################
def saveGen(generation, best):
    fname = os.path.join(gendir, "generation_%05d.pov" % generation)
    f = open(fname, 'w')
    f.write(repr(best))
    f.close()
    os.system("cd %s ; povray generation_%05d.pov > /dev/null 2>&1 &" % (gendir, generation))


###############################################################################
def doGeneration(gen, pop):
    global generation
    start = time.time()
    pop.gen()
    end = time.time()
    best = pop.best()
    saveGen(gen, best)
    bestfit = best.fitness()
    if bestfit <= 10:
        print "Finished at generation %d" % generation
        print best
        return True
    print "Generation %d (Bestfit %s) took %0.2fs" % (generation, bestfit, end-start)
    generation += 1


###############################################################################
def main():
    global generation
    try:
        os.mkdir(gendir)
    except OSError:
        pass
    pop = reflectorPopulation()
    reference_image.save("ref_img_scaled.png")

    try:
        for gen in range(1000):
            if doGeneration(gen, pop):
                break
        else:
            print "Never going to finish"
            print pop.best()
    except KeyboardInterrupt:
        print "Interrupted at generation %d" % generation

###############################################################################
if __name__ == "__main__":
    main()

#EOF
