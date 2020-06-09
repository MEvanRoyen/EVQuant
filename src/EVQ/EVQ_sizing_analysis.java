package EVQ;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.WaitForUserDialog;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;

/*
EVQuant FIJI Plugin
written by J.A. Slotman @ www.ErasmusOIC.nl
Based on FIJI macro's and excel templates By T.A. Hartjes
*/

public class EVQ_sizing_analysis implements PlugIn {

    private NumberFormat formatter = new DecimalFormat("#.###");

    public void run(String args){

        int Threshold = (int) Prefs.get("EVQ.Threshold",0);
        double SizeFactor = Prefs.get("EVQ.SizeCalibration",1);
        int SampleNum = (int) Prefs.get("EVQ.SizeSampleNum",1);
        int SampleRep = (int) Prefs.get("EVQ.SizeSampleRep",1);
        String[] model = {"Vesicle","Bead"};
        String modelChoice;


        GenericDialog gd = new GenericDialog("EVQuant Size analysis");
        gd.addNumericField("Threshold",Threshold,0);
        gd.addNumericField("Size Calibration Factor",SizeFactor,3);
        gd.addNumericField("Number of samples", SampleNum,0);
        gd.addNumericField("Number of stacks to analyse per sample", SampleRep,0);
        gd.addRadioButtonGroup("Used calibration model",model,1,2,model[0]);
        gd.showDialog();
        if(gd.wasCanceled()){
            return;
        }
        Threshold = (int) gd.getNextNumber();
        SizeFactor = gd.getNextNumber();
        SampleNum = (int) gd.getNextNumber();
        SampleRep = (int) gd.getNextNumber();
        modelChoice = gd.getNextRadioButton();

        getSizeCal sc = null;

        ResultsTable rt = new ResultsTable();


        for(int sample=0;sample<SampleNum;sample++){
            for(int rep=0;rep<SampleRep;rep++){

                WaitForUserDialog wt = new WaitForUserDialog("Open stack number "+(rep+1)+" from sample "+(sample+1));
                wt.show();
                ImagePlus imp = WindowManager.getCurrentImage();
                getEOST e = new getEOST(imp,Threshold,false);

                if(rep==0){
                    sc = new getSizeCal(imp,e,1);
                }else{
                    sc.addImp(e,imp);
                }

                WindowManager.getCurrentWindow().close();



            }

            sc.calculateSizes(SizeFactor,modelChoice);

            double meanSize = 0;


            for(Double i:sc.sizes){

                rt.incrementCounter();
                rt.addValue("Sample","Sample_"+(sample+1));
                rt.addValue("Particle size (diameter in nm)", i);
                meanSize += i;
            }

            meanSize /= sc.sizes.size();

            Plot p = new Plot("Size Distribution Sample_"+(sample+1),"Size of particles","Frequency");
            p.addHistogram(sc.sizesToArray(),4);
            p.show();
            IJ.log("Sample_"+(sample+1)+" Average Size : "+formatter.format(meanSize)+" Median Size : "+formatter.format(getMedian(sc.sizesToArray())));


        }

        rt.show("Results");

        Prefs.set("EVQ.Threshold",Threshold);
        Prefs.set("EVQ.SizeCalibration",SizeFactor);
        Prefs.set("EVQ.SizeSampleNum",SampleNum);
        Prefs.set("EVQ.SizeSampleRep",SampleRep);




    }

    public static double getMedian(double[] v){

        double[] values = v;
        double output;

        Arrays.sort(values);
        if (values.length % 2 != 0) {
            output = values[(int) Math.floor(values.length / 2.0)];
        } else {
            output = (  values[(int) Math.floor(values.length / 2.0)-1] + values[(int) Math.floor(values.length / 2.0)]   ) / 2.0;
        }

        return output;

    }

}
