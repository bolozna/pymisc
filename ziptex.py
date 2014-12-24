

helpstr="A tool to help you to make a zip package out of a .tex file and all the figures and the bibliography related to it.\nUsage: ziptex [texfile] \nOutputs a .zip file named such that foobar.tex produces foobar.zip."

import sys,os


if len(sys.argv)!=2:
    print helpstr
    exit(1)

inputfilename=sys.argv[1]
#inputfile=open(inputfilename,"r")

outputfilename=inputfilename[:-3]+"zip"

#figs=[]
#bibfile=None

def parseTex(inputfilename,figs,bibfiles,texfiles):
    inputfile=open(inputfilename,"r")
    for line in inputfile:
        line=line.strip()
        if not line.startswith("%"): #not a comment
            i=line.find(r"\includegraphics")
            if i!=-1:
                j=line[i:].find("{")
                k=line[j:].find("}")
                figs.add(line[i+j+1:j+k])

            i=line.find(r"\bibliography")
            if i!=-1:
                j=line[i:].find("{")
                k=line[j:].find("}")
                bibfiles.add(line[i+j+1:j+k])        

            i=line.find(r"\input{")
            if i!=-1:
                j=line[i:].find("{")
                k=line[j:].find("}")
                texfiles.add(line[i+j+1:j+k])        

unparsedFiles=set([inputfilename[:-4]])
figs=set()
bibfiles=set()
texfiles=[]
while len(unparsedFiles)>0:
    filename=unparsedFiles.pop()+".tex"
    print "Reading file: "+filename
    texfiles.append(filename)
    parseTex(filename,figs,bibfiles,unparsedFiles)


fileStr = " ".join(list(figs)) + " "+" ".join(list(texfiles))

print "In tex: Found "+str(len(figs))+" figures."
if len(bibfiles)>0:
    bibfileStr=" ".join(map(lambda x:x+".bib",list(bibfiles)))
    print "In tex: Found bibiliography: "+bibfileStr
    fileStr+=" "+bibfileStr

 

print "zip "+outputfilename+" "+fileStr
for line in os.popen("zip "+outputfilename+" "+fileStr): print line,

