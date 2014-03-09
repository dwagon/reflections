#!/usr/bin/env python
#
# Script to play around with genetic algorithms
#
import os
import sys
import getopt
import tempfile
import time
from pygene.population import Population
from pygene.organism import MendelOrganism
from pygene.gene import FloatGene
import Image
import ImageChops
import ImageStat

opts = {
    'single': False,
    }

verbFlag = False
worldWidth = 20
worldHeight = 20
worldDepth = 20
imageWidth = 64
imageHeight = 64
gendir = "gendir-%d" % os.getpid()
numreflectors = 20
stats = {'analyse': 0.0, 'generate': 0.0, 'render': 0.0}

ri = Image.open("reference.png")
reference_image = ri.resize((imageWidth, imageHeight))

generation = 0
renders = 0


###############################################################################
def Verbose(msg):
    if verbFlag:
        sys.stderr.write("%s\n" % msg)


###############################################################################
def Warning(msg):
    sys.stderr.write("Warning: %s\n" % msg)


###############################################################################
def Fatal(msg):
    sys.stderr.write("Fatal: %s\n" % msg)
    sys.exit(255)


###############################################################################
def usage():
    sys.stderr.write("Usage: %s\n" % sys.argv[0])


###############################################################################
class colourGene(FloatGene):
    mutProb = 0.05
    randMin = 0.0
    randMax = 1.0

    ###########################################################################
    def __repr__(self):
        return "C%0.2f" % self.value


###############################################################################
class heightGene(FloatGene):
    mutProb = 0.2
    randMin = 0.0
    randMax = worldHeight

    ###########################################################################
    def __repr__(self):
        return "H%0.2f" % self.value


###############################################################################
class depthGene(FloatGene):
    mutProb = 0.2
    randMin = 0.0
    randMax = worldDepth

    ###########################################################################
    def __repr__(self):
        return "D%0.2f" % self.value


###############################################################################
class radiusGene(FloatGene):
    mutProb = 0.1
    randMin = 0.5
    randMax = 1.0

    ###########################################################################
    def __repr__(self):
        return "R%0.2f" % self.value


###############################################################################
class widthGene(FloatGene):
    mutProb = 0.2
    randMin = 0.0
    randMax = worldWidth

    ###########################################################################
    def __repr__(self):
        return "W%f" % self.value

genome = {}
for i in range(numreflectors):
    genome["x_%d" % i] = widthGene
    genome["y_%d" % i] = heightGene
    genome["z_%d" % i] = depthGene
    genome["r_%d" % i] = radiusGene
    genome["cr_%d" % i] = colourGene
    genome["cg_%d" % i] = colourGene
    genome["cb_%d" % i] = colourGene


###############################################################################
class reflectBox(MendelOrganism):
    genome = genome

    ###########################################################################
    def __repr__(self):
        delta = 0.01
        str = ""
        str += '#include "colors.inc"\n'
        str += '#include "finish.inc"\n'
        str += 'camera { location <-0.5, -0.5, -0.5> look_at  <9, 9, 9> angle 90 }\n'
        str += 'global_settings { ambient_light White }\n'
        str += 'light_source { <%s, %s, %s> color White }\n' % (-1+delta, -1+delta, -1+delta)
        str += 'light_source { <%s, %s, %s> color White }\n' % (worldWidth-delta, -1+delta, -1+delta)
        str += 'light_source { <%s, %s, %s> color White }\n' % (worldWidth-delta, -1+delta, worldDepth+delta)
        str += 'light_source { <%s, %s, %s> color White }\n' % (worldWidth-delta, worldHeight-delta, worldDepth-delta)
        # str += 'plane { <0, 1, 0>, %s pigment { Black }  }\n' % (worldHeight+1)    # top
        str += 'plane { <0, 1, 0>, -1 pigment { White }  }\n'   # bottom
        # str += 'plane { <1, 0, 0>, -1 pigment { Black }  }\n'  # left
        # str += 'plane { <1, 0, 0>, %s pigment { White }  }\n' % (worldWidth+1)   # right
        # str += 'plane { <0, 0, 1>, -1 pigment { Black }  }\n'  # front
        str += 'plane { <0, 0, 1>, %s pigment { Yellow }  }\n' % (worldDepth+1)   # back
        for o in range(numreflectors):
            colour = "rgb <%f, %f, %f>" % (self["cr_%d" % o], self["cg_%d" % o], self["cb_%d" % o])
            sphere = "sphere { <%f, %f, %f>, %f texture { pigment { color %s } }}\n" % (self["x_%d" % o], self["y_%d" % o], self["z_%d" % o], self["r_%d" % o], colour)
            str += sphere
        return str

    ###########################################################################
    def generateScene(self):
        start = time.time()
        f = tempfile.NamedTemporaryFile(prefix="%s-" % generation, delete=False, suffix='.pov')
        f.write(str(self))
        f.close()
        stats['generate'] += (time.time()-start)
        return f.name

    ###########################################################################
    def renderScene(self, fname):
        start = time.time()
        pngf = fname.replace('.pov', '.png')
        f = os.popen('povray -D +H%s +W%s +O%s %s 2>&1' % (imageHeight, imageWidth, pngf, fname))
        log = f.read()
        x = f.close()
        if x:
            sys.stderr.write(log)
        stats['render'] += (time.time()-start)
        return pngf

    ###########################################################################
    def analyseScene(self, fname):
        # Lower fitness is better
        start = time.time()
        image = Image.open(fname)
        diffimage = ImageChops.difference(reference_image, image)
        diff = sum([x for x in ImageStat.Stat(diffimage).rms])
        diff = diff*diff*diff   # Make differences big
        stats['analyse'] += (time.time()-start)
        return diff

    ###########################################################################
    def fitness(self):
        global renders
        renders += 1
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
class reflectorPopulation(Population):
    initPopulation = 20
    species = reflectBox
    childCull = 2
    childCount = 10
    mutants = 0.20
    numNewOrganisms = 2
    incest = 1


###############################################################################
def saveGen(generation, best):
    fname = os.path.join(gendir, "generation_%05d.pov" % generation)
    f = open(fname, 'w')
    f.write(repr(best))
    f.close()
    os.system("cd %s ; povray generation_%05d.pov > /dev/null 2>&1 &" % (gendir, generation))


###############################################################################
def doGeneration(gen, pop):
    global renders, stats, generation
    genfitness = 0
    genlength = 0
    start = time.time()
    b = pop.best()
    saveGen(gen, b)
    bestfit = b.fitness()
    if bestfit <= 1:
        print "Finished at generation %d" % generation
        print b
        return
    pop.gen()
    print "Generation %d x %d (%s) took %0.2fs [generate %0.2fs render %0.2fs (%d) analyse %0.2fs]" % (generation, genlength, bestfit, time.time()-start, stats['generate'], stats['render'], renders, stats['analyse'])
    renders = 0
    if bestfit != genfitness:
        genfitness = bestfit
        genlength = 1
    else:
        genlength += 1
    if generation > 100 and genlength > 30:
        print "Not moving anywhere"
        return
    generation += 1
    stats = {'analyse': 0.0, 'generate': 0.0, 'render': 0.0}


###############################################################################
def main():
    global generation, stats
    try:
        os.mkdir(gendir)
    except OSError:
        pass
    stats = {'analyse': 0.0, 'generate': 0.0, 'render': 0.0}
    pop = reflectorPopulation()
    reference_image.save("ref_img_scaled.png")

    try:
        for gen in range(1000):
            doGeneration(gen, pop)
    except KeyboardInterrupt:
        print "Interrupted at generation %d" % generation
    else:
        print "Never going to finish"
        b = pop.best()
        print b

###############################################################################
if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "vh", ["single"])
    except getopt.GetoptError, err:
        sys.stderr.write("Error: %s\n" % str(err))
        usage()
        sys.exit(1)

    for o, a in opts:
        if o == "-v":
            verbFlag = True
        if o == "-h":
            usage()
            sys.exit(0)
        if o == "--single":
            opts['single'] = True

    main()

#EOF
