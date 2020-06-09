package EVQ;


import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.measure.CurveFitter;
import ij.measure.ResultsTable;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/*
EVQuant FIJI Plugin
written by J.A. Slotman @ www.ErasmusOIC.nl
Based on FIJI macro's and excel templates By T.A. Hartjes
*/

public class getSizeCal {

    ArrayList<Double> hGram = new ArrayList<Double>();
    ArrayList<Double> sizes = new ArrayList<Double>();
    int BinWidth;
    double beadSize;
    double volcor;
    double volcor_bead;
    double[] volcor_results;
    static final int HALFLINE = 8;
    static final int LINESIZE = (HALFLINE*2)+1;
    static final double FWHMFACTOR = 2 * Math.sqrt(2*Math.log(2));

    public getSizeCal(ImagePlus imp, getEOST e, double beadSize_){

        beadSize = beadSize_;

        getVolumeHisto(e,imp);
        measureVolcor(beadSize);

    }


    public void showResults(){
        ResultsTable rt = new ResultsTable();
        for(Double i : hGram){
            rt.incrementCounter();
            rt.addValue("Intensity",i);
        }

        rt.show("Raw Results Sizing analysis");
    }


    public void addImp(getEOST e, ImagePlus imp){
        getVolumeHisto(e,imp);
        measureVolcor(beadSize);

    }

    public void measureVolcor(double beadSize){

        double meanSize = 0;

        for(Double i:hGram){
            meanSize += i;
        }

        meanSize /= hGram.size();



        Plot p = new Plot("Histogram of particle sizes","Intensity (Binsize: "+BinWidth+")","Frequency");
        p.addHistogram(hGramToArray(),BinWidth);

        float[] xf = p.getXValues();
        float[] yf = p.getYValues();

        double[] x = new double[xf.length];
        double[] y = new double[yf.length];

        for(int i=0;i<x.length;i++){
            x[i] = (double) xf[i];
            y[i] = (double) yf[i];
        }

        CurveFitter cv = new CurveFitter(x,y);
        cv.doFit(cv.GAUSSIAN);
        double[] results = cv.getParams();
        volcor_results = results;

        volcor = results[2] / Math.pow(beadSize/2,2);
        volcor_bead = results[2] / Math.pow(beadSize/2,3);

    }

    public void calculateSizes(double calib, String model){
        if(model.equals("Vesicle")) {
            for (Double i : hGram) {
                double size = Math.sqrt(i / calib) * 2;
                sizes.add(size);
            }
        }
        if(model.equals("Bead")){
            for (Double i : hGram) {
                double size = Math.cbrt(i / calib) * 2;
                sizes.add(size);
            }
        }
    }

    private double[] hGramToArray(){

        double[] out = new double[hGram.size()];

        for(int i=0;i<hGram.size();i++){
            out[i] = hGram.get(i);
        }

        return out;
    }

    public double[] sizesToArray(){

        double[] out = new double[sizes.size()];

        for(int i=0;i<sizes.size();i++){
            out[i] = sizes.get(i);
        }

        return out;
    }

    public void plot(){

        Plot p = new Plot("Histgram of particle sizes","Intensity (Binsize: "+BinWidth+")","Frequency");
        p.addHistogram(hGramToArray(),BinWidth);
        float[] xvalues_raw = p.getXValues();
        double[] xvalues = new double[xvalues_raw.length];
        double[] yvalues = new double[xvalues_raw.length];
        for(int i=0;i<xvalues.length;i++){

            xvalues[i] = (double) xvalues_raw[i];
            yvalues[i] =  volcor_results[0] + (volcor_results[1]-volcor_results[0])*Math.exp( - Math.pow(xvalues[i]-volcor_results[2],2)/(2*Math.pow(volcor_results[3],2)) );

        }
        p.setColor(Color.RED);
        p.add("line",xvalues,yvalues);
        p.show();

    }


    public void getVolumeHisto(getEOST e, ImagePlus imp){

        Calibration c = imp.getCalibration();
        int[][] centerPoints = e.centerPoints;
        ImageStack ims = imp.getImageStack();
        ResultsTable rt = new ResultsTable();
        boolean checkFit = true;


        for(int[] point: centerPoints) {

            IJ.showStatus("Performing Size Calibration");


            if (point[2] > HALFLINE  && point[2] < ims.size() - HALFLINE) {

                if(checkFit) {
                    rt.incrementCounter();
                }else{
                    checkFit = true;
                }

                double[][] index = new double[3][LINESIZE];
                double[] xProfile = new double[LINESIZE];
                double[] yProfile = new double[LINESIZE];
                double[] zProfile = new double[LINESIZE];



                ImageProcessor ip = ims.getProcessor(point[2]);
                ImageProcessor ipz;

                for (int i = 0; i < LINESIZE; i++) {
                    index[0][i] = point[0] - HALFLINE + i;
                    index[1][i] = point[1] - HALFLINE + i;
                    index[2][i] = point[2] - HALFLINE + i;
                    xProfile[i] = ip.getPixel((int)index[0][i], point[1]);
                    yProfile[i] = ip.getPixel(point[0], (int)index[1][i]);
                    ipz = ims.getProcessor((int) index[2][i] );
                    zProfile[i] = ipz.getPixel(point[0], point[1]);
                }

                double fwhmZ = 0;
                double[] ind = null;
                double[] yvalues = new double[LINESIZE];
                int size = 0;

                for (int i = 0; i < 3; i++) {

                    switch (i) {
                        case 0:
                            yvalues = xProfile;
                            size = ims.getWidth();
                            ind = index[0];
                            break;
                        case 1:
                            yvalues = yProfile;
                            size = ims.getHeight();
                            ind = index[1];
                            break;
                        case 2:
                            yvalues = zProfile;
                            size = ims.getSize();
                            ind = index[2];
                            break;
                    }

                    CurveFitter cv = new CurveFitter(ind, yvalues);
                    cv.doFit(CurveFitter.GAUSSIAN);
                    double[] results = cv.getParams();
                    double goodness = cv.getRSquared();

                    for (int j = 0; j < 5; j++) {
                        rt.addValue("Param_" + i + "_" + j, results[j]);
                    }

                    rt.addValue("Param_"+i+"_5",goodness);

                    if(results[2]<0 || results[2]>size  || results[1]>Math.pow(2,16)){
                        checkFit = false;
                    }



                }





            }
        }

        //remove last row if it was outside fit parameters
        if(!checkFit){
            rt.deleteRow(rt.getCounter()-1);
        }



        ResultsTable rt2 = new ResultsTable();

        double meanSize = 0;
        double z = 0;
        double maxz = 0;

        double background = imp.getStatistics().mean;

        for(int i=0;i<rt.size();i++){

            z = rt.getValue("Param_2_2",i);
            meanSize = (rt.getValue("Param_0_1",i) + rt.getValue("Param_1_1",i) + rt.getValue("Param_2_1",i))/3 ;
            maxz = rt.getValue("Param_2_1",i);

            rt2.incrementCounter();
            rt2.addValue("z",z);
            rt2.addValue("mean Max (xyz)",meanSize);
            rt2.addValue("Max z",maxz-background);
        }


        //rt2.show("Results");

        double maxHgram = 0;
        double minHgram = Math.pow(2,16);

        for(int i=0;i<rt.size();i++){

            double v = rt2.getValue("Max z",i);

            if(v < Math.pow(2,16) && v > 0) {
                hGram.add(v);
                if(v>maxHgram){
                    maxHgram = v;
                }
                if(v<minHgram){
                    minHgram = v;
                }
            }

        }

        BinWidth = (int) getFDbinsize( hGramToArray() );



    }


    public static double getFDbinsize(double[] v){

        double[] values = v;
        double Q1;
        double Q2;
        double Q3;

        Arrays.sort(values);
        double q4 = values.length;
        double q1,q2,q3;

        q2 = q4/2.0;

        if(q2 % 2 != 0){
            q1 = Math.floor(q2) / 2.0;
            q3 = q4 - q1;
            Q2 = values[ (int) Math.floor(q2) ];
        }else{
            q1 = q2 / 2.0;
            q3 = q4 - q1;
            Q2 = (values[(int)(q1-1)] + values[(int)q1])/2.0;
        }

        if(q1 % 2 != 0){
            Q1 = values[ (int) Math.floor(q1)];
        }else{
            Q1 = values[(int)(q1-1)] + values[(int)q1];
        }

        if(q3 % 2 != 0){
            Q3 = values[ (int) Math.floor(q3)];
        }else{
            Q3 = values[(int)(q3-1)] + values[(int)q3];
        }


        double binsize = 2 * ((Q3-Q1)/ Math.cbrt(q4));

        return Math.floor(binsize) + 1;
    }



}
