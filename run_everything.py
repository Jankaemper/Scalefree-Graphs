import subprocess
import sys
import argparse

def compile():
    subprocess.Popen(
    "cc -o shortest_path shortest_path_fragment.c shortest_path_main.c list.c -lm",
    shell=True).wait()

def simulationScalefree():
    subprocess.Popen("./ssimulationScalefree", shell=True).wait()

def simulationPlanar():
    subprocess.Popen("./simulationPlanar.scr", shell=True).wait()

def simulationConstant():
    subprocess.Popen("./simulationConstant.scr", shell=True).wait()

def manual_run():
    print("mode: 1-Scale Free 2-Scale Free Constant Pref 3-Planar 4-Planar Mean\n");
    mode = int(input("Welcher mode soll gewaehlt werden?: "))
    if mode == 1:
        node = int(input("Anzahl der Knoten: "))
        m  =  int(input("M: "))
        eingabe = "./shortest_path " + str(mode) +" "+str(node)+ " " + str(m)
    if mode == 2:
        node = int(input("Anzahl der Knoten: "))
        m  =  int(input("M: "))
        print("Bedingung: k_0 > -m")
        print("k_0 als double angeben")
        k_0 = (input("k_0: "))
        eingabe = "./shortest_path " + str(mode) +" "+str(node)+ " " + str(m) +" "+ str(k_0)
    if mode == 3:
        node = int(input("Anzahl der Knoten: "))
        eingabe = "./shortest_path " + str(mode) +" "+str(node)
    if mode == 4:
        maximum = int(input("max: "))
        stepSize =  int(input("stepSize: "))
        eingabe = "./shortest_path " + str(mode) +" "+str(maximum)+ " " + str(stepSize)
    print eingabe
    print ("Calculating ...")
    subprocess.Popen(eingabe, shell=True).wait()

def gnuplotConstant():
    subprocess.Popen("gnuplot Plotscripts/gnuplot_ausgabe_fit_constant_powerLaw_M2.plt -", shell=True).wait()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--compile", help="Compile into shortest_path", action="store_true")
    parser.add_argument("-simS", "--simulationScalefree", help="run ./simulationScalefree.src", action="store_true")
    parser.add_argument("-simC", "--simulationConstant", help="run ./simulationConstant.src", action="store_true")
    parser.add_argument("-simP", "--simulationPlanar", help="run ./simulationPlanar.src", action="store_true")
    parser.add_argument("-m", "--manual", help="chose parameter manually", action="store_true")
    parser.add_argument("-gnuC", "--gnuplotConstant", help="create diagramms for Constant", action="store_true")
    args = parser.parse_args()
    if len(sys.argv) == 1:
        print("** shortest_path calculation  **")
        print("** For USAGE use run_everything.py -h **")
    if args.manual:
        manual_run()
    if args.compile:
        compile()
    if args.simulationConstant:
        simulationConstant()
    if args.simulationScalefree:
        simulationScalefree()
    if args.simulationPlanar:
        simulationPlanar()
    if args.gnuplotConstant:
        gnuplotConstant()
