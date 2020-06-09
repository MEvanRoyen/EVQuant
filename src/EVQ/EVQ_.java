package EVQ;

import ij.*;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.WaitForUserDialog;
import ij.plugin.PlugIn;
import ij.plugin.filter.GaussianBlur;
import ij.process.ImageProcessor;


import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;

/*
EVQuant FIJI Plugin
written by J.A. Slotman @ www.ErasmusOIC.nl
Based on FIJI macro's and excel templates By T.A. Hartjes

Version 1.0.2096 : Initial upload to github
*/

public class EVQ_ implements PlugIn {


    private static final double detectionThreshold = 1.0/2500.0; // 1 particle per 2500 um^2
    int threshold;
    private NumberFormat formatter = new DecimalFormat("#.###");
    private NumberFormat formatter2 = new DecimalFormat("#.######");

    public void run(String args) {

        boolean doSizing = false;


        if(args.equals("all")) {

            GenericDialog startMenu = new GenericDialog("EVQuant Calibration");
            startMenu.addMessage("Threshold calibration");
            startMenu.addMessage("Estimated Optical Slice Thickness (EOST) calibration");
            startMenu.addMessage("______________________________________________________");
            startMenu.addMessage("Optional calibrations:");
            startMenu.addCheckbox("Particle sizing calibration", doSizing);


            startMenu.showDialog();

            if (startMenu.wasCanceled()) {
                return;
            }

            doSizing = startMenu.getNextBoolean();


        }

        if(args.equals("calibration") || args.equals("all")) {

            GenericDialog gd = new GenericDialog("EVQuant Calibration");
            gd.addNumericField("Number of samples for calibration curve including negative control", 4, 0, 3, null);
            gd.addCheckbox("Define alternative threshold as a starting point",false);
            gd.showDialog();

            if (gd.wasCanceled()) {
                return;
            }

            int CalibNumber = (int) gd.getNextNumber();
            boolean altstart = gd.getNextBoolean();
            int startingPoint = 0;




            double[] concentrations = new double[CalibNumber];
            double[] sds = new double[CalibNumber];
            MyOwnFindmaxima mf = new MyOwnFindmaxima();

            int outputType = 5;

            double[][] count = new double[CalibNumber][255];


            for (int c = 0; c < CalibNumber; c++) {

                WaitForUserDialog wt;
                if (c == 0) {
                    wt = new WaitForUserDialog("Please open the stack of images of the negative control sample and select OK when the file is loaded");
                } else {
                    wt = new WaitForUserDialog("Please open the stack of images of dilution number " + c +" and select OK when the file is loaded");
                }
                wt.show();


                ImagePlus imp = WindowManager.getCurrentImage();

                if(imp==null){
                    IJ.showMessage("error","No image is open");
                    return;
                }

                if(altstart && c==0){
                    GenericDialog gd_altstart = new GenericDialog("EVQuant Calibration");
                    gd_altstart.addNumericField("Alternative threshold used as a start point",0,0);
                    gd_altstart.showDialog();
                    if (gd_altstart.wasCanceled()) {
                        return;
                    }

                    startingPoint = (int) gd_altstart.getNextNumber();

                    if(startingPoint>= Math.pow(2,imp.getBitDepth())){
                        IJ.showMessage("error","Threshold is bigger than highest possible intensity");
                        return;
                    }
                    if(startingPoint<0){
                        IJ.showMessage("error","Threshold is smaller than zero");
                        return;
                    }


                }

                if(c==0){
                    count = new double[CalibNumber][(int)Math.pow(2,imp.getBitDepth())];
                }

                double imageArea = (imp.getWidth()*imp.getCalibration().pixelWidth) * (imp.getHeight()*imp.getCalibration().pixelHeight);

                ImageStack ims = imp.getStack();

                if(ims.size()<=1){
                    IJ.showMessage("error", "Image is not a stack");
                    return;
                }

                if (c != 0) {
                    GenericDialog gd2 = new GenericDialog("EVQuant Calibration");
                    gd2.addNumericField("Concentration dilution number " + c , 0, 1, 3, null);
                    gd2.showDialog();
                    if (gd2.wasCanceled()) {
                        return;
                    }
                    concentrations[c] = gd2.getNextNumber();
                } else {
                    concentrations[c] = 0;
                }


                ImageProcessor[] ips = new ImageProcessor[ims.size()];
                GaussianBlur gb = new GaussianBlur();

                for (int stack = 1; stack <= ims.size(); stack++) {
                    ips[stack - 1] = ims.getProcessor(stack);
                    gb.blurGaussian(ips[stack - 1], 2, 2, 0.002);
                }

                if (c == 0) {

                    for (int i = startingPoint; i < count[0].length; i++) {

                        int[] values = new int[ims.size()];


                        for (int stack = 0; stack < ims.size(); stack++) {
                            values[stack] = (int) mf.findMaxima(ips[stack], i, outputType, true);

                        }

                        sds[c] = getSD(values);
                        count[c][i] = getMedian(values);

                        IJ.showStatus("Calculating threshold "+i);

                        if (count[c][i] <= imageArea * detectionThreshold || count[c][i] <= 1 ) {
                            threshold = i;
                            i = count[0].length;
                        }



                    }


                } else {

                    int[] values = new int[ims.size()];

                    for (int stack = 0; stack < ims.size(); stack++) {
                        values[stack] = (int) mf.findMaxima(ips[stack], threshold, outputType, true);
                    }

                    sds[c] = getSD(values);
                    count[c][threshold] = getMedian(values);
                }

                WindowManager.closeAllWindows();




            }



            IJ.log("Calibrated threshold " + threshold);

            double[] yvalues = new double[concentrations.length];
            double xmin = concentrations[0];
            double xmax = concentrations[0];
            double ymin = yvalues[0];
            double ymax = yvalues[0];


            for (int i = 0; i < CalibNumber; i++) {
                yvalues[i] = count[i][threshold];
                IJ.log(concentrations[i] + " , " + count[i][threshold]);

                if (yvalues[i] - sds[i] < ymin) {
                    ymin = yvalues[i] - sds[i];
                }
                if (yvalues[i] + sds[i] > ymax) {
                    ymax = yvalues[i] + sds[i];
                }
                if (concentrations[i] < xmin) {
                    xmin = concentrations[i];
                }
                if (concentrations[i] > xmax) {
                    xmax = concentrations[i];
                }

            }

            for(int i=0;i<concentrations.length-1;i++){
                for(int j=i+1;j<concentrations.length;j++){

                    if(concentrations[i]>concentrations[j]){
                        double swap = concentrations[i];
                        concentrations[i] = concentrations[j];
                        concentrations[j] = swap;

                        swap = yvalues[i];
                        yvalues[i] = yvalues[j];
                        yvalues[j] = swap;

                    }


                }
            }



            Plot plot = new Plot("EVQuant Calibration Curve", "Concentration", "Vesicle Counts");
            plot.setLimits(xmin - ((xmax - xmin) / 10.0), xmax + ((xmax - xmin) / 10.0), ymin - ((ymax - ymin) / 10.0), ymax + ((ymax - ymin) / 10.0));
            plot.add("line", concentrations, yvalues);
            plot.addPoints(concentrations, yvalues, sds, Plot.CIRCLE);
            plot.addLegend("Error Bars\tStandard Deviation");

            plot.show();

            Prefs.set("EVQ.Threshold",threshold);
            Prefs.savePreferences();

            WindowManager.setWindow( WindowManager.getWindow("Log") );


        }

        if(args.equals("eost") || args.equals("all")){

            int thr = (int) Prefs.get("EVQ.Threshold",0);
            int stacknumber = 1;

            boolean report = true;

            GenericDialog gd_eost = new GenericDialog("Effictive Optical Slice Thickness (EOST) calibration");
            gd_eost.addNumericField("Number of stacks used for EOST calibration",stacknumber,0);


            if(args.equals("all")) {

                thr = threshold;

                report = false;


            }

            if(!args.equals("all")) {
                gd_eost.addNumericField("Enter Threshold", thr, 0);
            }

            gd_eost.showDialog();

            if (gd_eost.wasCanceled()) {
                   return;
            }

            stacknumber = (int) gd_eost.getNextNumber();

            if(!args.equals("all")) {
                thr = (int) gd_eost.getNextNumber();
            }



            double eosts =0;
            double[] eostsArray = new double[stacknumber];


            for(int i=0;i<stacknumber;i++) {

                WaitForUserDialog wt = new WaitForUserDialog("Please open calibration stack number "+(i+1)+" for EOST measurement");
                wt.show();

                ImagePlus imp = WindowManager.getCurrentImage();
                getEOST e = new getEOST(imp, thr, report);


                if(args.equals("all")) {
                    WindowManager.getCurrentWindow().close();
                }

                eostsArray[i] = e.EOST;
                eosts += e.EOST;






            }

            eosts = eosts/stacknumber;

            double sd = 0;
            if(stacknumber>1) {
                sd = getSD(eostsArray);
            }



            IJ.log("Average EOST: " + formatter.format(eosts) +" sd: "+formatter.format(sd));

            Prefs.set("EVQ.EOST",eosts);
            //Prefs.savePreferences();

            IJ.showProgress(1);
            IJ.showStatus("Finished...");

            WindowManager.setWindow( WindowManager.getWindow("Log") );




        }

        if(args.equals("Sizing") || doSizing){


            int thr = (int) Prefs.get("EVQ.Threshold",0.0);
            boolean report = false;
            double beadSize = Prefs.get("EVQ.SizeBeadsize",50);
            int repnum = (int) Prefs.get("EVQ.SizeRepnum",1);
            String[] model = {"Vesicle","Bead"};
            String modelChoice;
            boolean results = false;
            boolean plot = false;

            GenericDialog gd_sizing = new GenericDialog("Particle sizing calibration");
            if(!args.equals("all")){
                gd_sizing.addNumericField("Threshold",thr,0);
            }
            gd_sizing.addNumericField("Known particle size of calibration sample (diameter in nm)",beadSize,1);
            gd_sizing.addNumericField("Number of stacks used for sizing calibration",repnum,0);
            gd_sizing.addRadioButtonGroup("Calibration Model",model,1,2,model[0]);
            gd_sizing.addCheckbox("Show Raw Results",results);
            gd_sizing.addCheckbox("Show plot",plot);

            gd_sizing.showDialog();

            if(gd_sizing.wasCanceled()){
                return;
            }

            if(!args.equals("all")){
                thr = (int) gd_sizing.getNextNumber();
            }

            beadSize = gd_sizing.getNextNumber();
            repnum = (int) gd_sizing.getNextNumber();
            modelChoice = gd_sizing.getNextRadioButton();
            results = gd_sizing.getNextBoolean();
            plot = gd_sizing.getNextBoolean();
            getSizeCal sc = null;

            for(int i=0;i<repnum;i++) {

                WaitForUserDialog wt = new WaitForUserDialog("Please open calibration stack number "+(i+1)+" for sizing analysis");
                wt.show();
                ImagePlus imp = WindowManager.getCurrentImage();

                IJ.showStatus("Determining centers of particles");
                getEOST e = new getEOST(imp, thr, report);
                if(i==0) {
                    sc = new getSizeCal(imp, e, beadSize);
                }else{
                    sc.addImp(e,imp);
                }

                WindowManager.closeAllWindows();

            }

            if(plot) {
                sc.plot();
            }
            if(results) {
                sc.showResults();
            }

            double factor=0;

            if(modelChoice.equals("Vesicle")){
                factor = sc.volcor;
            }
            if(modelChoice.equals("Bead")){
                factor = sc.volcor_bead;
            }

            IJ.log("Size Calibration Factor: "+formatter2.format(factor));


            Prefs.set("EVQ.SizeRepnum",repnum);
            Prefs.set("EVQ.SizeBeadsize",beadSize);
            Prefs.set("EVQ.SizeCalibration",factor);

            WindowManager.setWindow( WindowManager.getWindow("Log") );

        }

    }



    public double getSD(int[] v){
        double sd=0;
        double mean=0;

        for(double x : v){
            mean += x;
        }

        mean /= v.length;

        for(double x : v){
            sd += Math.pow(x-mean,2);
        }

        sd = Math.sqrt( sd/ (v.length - 1) );

        return sd;
    }

    public double getSD(double[] v){
        double sd=0;
        double mean=0;

        for(double x : v){
            mean += x;
        }

        mean /= v.length;

        for(double x : v){
            sd += Math.pow(x-mean,2);
        }

        sd = Math.sqrt( sd/ (v.length - 1) );

        return sd;
    }

    public double getMedian(int[] v){

        int[] values = v;
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
