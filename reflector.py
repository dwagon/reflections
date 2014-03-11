#!/usr/bin/env python
#
# Script to play around with genetic algorithms
#
import os
import sys
import tempfile
import time
import math
from itertools import izip
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
numreflectors = 1

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
    randMax = 10.0

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
        genome["r_%d" % i] = radiusGene
        #genome["cr_%d" % i] = colourGene
        #genome["cg_%d" % i] = colourGene
        #genome["cb_%d" % i] = colourGene

    ###########################################################################
    def __repr__(self):
        delta = 0.01
        wwd = worldWidth - delta
        wdd = worldDepth - delta
        str = ""
        str += '#include "colors.inc"\n'
        str += '#include "finish.inc"\n'
        str += '#include "textures.inc"\n'
        str += 'camera { location <-0.5, -0.5, -0.5> look_at  <9, 9, 9> angle 90 }\n'
        str += 'global_settings { ambient_light White }\n'
        #str += 'light_source { <%s, %s, %s> color White }\n' % (-1+delta, -1+delta, -1+delta)
        #str += 'light_source { <%s, %s, %s> color White }\n' % (wwd, -1+delta, -1+delta)
        #str += 'light_source { <%s, %s, %s> color White }\n' % (wwd, -1+delta, worldDepth+delta)
        #str += 'light_source { <%s, %s, %s> color White }\n' % (wwd, worldHeight-delta, wdd)
        #str += 'plane { <0, 1, 0>, %s pigment { Black }  }\n' % (worldHeight+1)    # top
        str += 'plane{ <1,0,0>, -2 texture{ pigment{ checker color rgb<1,1,1>*1.2 color rgb<0.25,0.15,0.1>*0} normal { bumps 0.75 scale 0.025} finish { phong 0.1}}}\n'
        str += 'plane{ <0,1,0>, -2 texture{ pigment{ checker color rgb<1,1,1>*1.2 color rgb<0.1,0.25,0.15>*0} normal { bumps 0.75 scale 0.025} finish { phong 0.1}}}\n'
        str += 'plane{ <0,0,1>, -2 texture{ pigment{ checker color rgb<1,1,1>*1.2 color rgb<0.1,0.15,0.25>*0} normal { bumps 0.75 scale 0.025} finish { phong 0.1}}}\n'
        #str += 'plane { <1, 0, 0>, %s finish { reflection { 0.75 } ambient 0.1 }}\n' % (worldDepth+1) 
        #str += 'plane { <1, 0, 0>, %s pigment { White }  }\n' % (worldWidth+1)   # right
        #str += 'plane { <0, 0, 1>, -1 pigment { Black }  }\n'  # front
        str += 'plane { <0, 0, 1>, %s finish { reflection { 0.75 } ambient 0.1 }}\n' % (worldDepth+1) 
        for o in range(numreflectors):
            colour = "rgb <9.0, 0.0, 0.0>"
            #colour = "rgb <%f, 0.0, 0.0>" % (self["cr_%d" % o], )
            sphere = "sphere { <9.0, 9.0, 9.0>, %f texture { Polished_Chrome pigment {color %s } }}\n" % (self["r_%d" % o], colour)
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
        pairs = izip(image.getdata(), reference_image.getdata())
        dif = sum(abs(c1-c2) for p1,p2 in pairs for c1,c2 in zip(p1,p2))
        radiussum=0
        for o in range(numreflectors):
            radiussum += 10.0*self["r_%d" % o]
        return dif + radiussum

    ###########################################################################
    def fitness(self):
        povfile = self.generateScene()
        pngfile = self.renderScene(povfile)
        ftn = self.analyseScene(pngfile)
        self.cleanup(povfile, pngfile)
        f=open('fitness.log','a')
        f.write("%f,%f\n" % (self['r_0'], ftn))
        f.close()
        return ftn

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
