package EVQ;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.measure.CurveFitter;
import ij.measure.ResultsTable;


import ij.process.ImageProcessor;

import java.awt.*;
import java.util.ArrayList;

/*
EVQuant FIJI Plugin
written by J.A. Slotman @ www.ErasmusOIC.nl
Based on FIJI macro's and excel templates By T.A. Hartjes
*/

public class EVQ_imunno  {

    int R = 3;
    int D = (2*R)+1;
    int BGR = 15;
    int BGD = (2*BGR)+1;
    int CIRCLEAREA = new OvalRoi(0,0, D, D).getContainedPoints().length;
    int BGCIRCLEAREA = new OvalRoi(0,0, BGD, BGD).getContainedPoints().length;
    boolean positives[][];
    double[] thrs;
    double[] medians;

    String[] CollumnNames;
    String[] ChNames;
    double[] Values;

    public void EVQ_immuno(){


    }

    public void immuno(ImagePlus[] impCh, Calibration c, int refChNo, int threshold, boolean report, String[] channelNames,
                       double eots,int totalVolume, int dups, double sampleVolume, double dil, double[] cutoffs,String sample, int edgesize ){

        ResultsTable rt = new ResultsTable();


        getEOST e = new getEOST(impCh[refChNo],threshold,false);
        int nCh = impCh.length;

        double[][] averages_r = new double[e.allPoints.length][nCh];
        double[][] averages_r_dil = new double[e.allPoints.length][nCh];



        int cnt = 0;
        for(int[] point : e.allPoints ){

            if(point[0]-(31+edgesize) > 0 && point[0]+(31+edgesize) < impCh[0].getWidth() && point[1]-(31+edgesize) > 0 && point[1]+(31+edgesize) < impCh[0].getHeight()) {

                OvalRoi r = new OvalRoi(point[0] - R, point[1] - R, D, D);
                OvalRoi r_dilated = new OvalRoi(point[0] - BGR, point[1] - BGR, BGD, BGD);

                Point[] roiPix = r.getContainedPoints();
                Point[] roiPix_dil = r_dilated.getContainedPoints();

                for (int i = 0; i < nCh; i++) {
                    ImageProcessor ip = impCh[i].getImageStack().getProcessor(point[2]);

                    for (Point p : roiPix) {
                        averages_r[cnt][i] += ip.get((int) p.getX(), (int) p.getY());
                    }

                    for (Point p : roiPix_dil) {

                        //IJ.log("1: "+p.getX()+"2: "+p.getY()+" roi "+cnt );
                        averages_r_dil[cnt][i] += ip.get((int) p.getX(), (int) p.getY());
                    }




                }

                cnt++;
            }


        }

        for(int i=0;i<cnt;i++){

            rt.incrementCounter();

            for(int j=0; j<averages_r[0].length;j++) {

                rt.addValue("Area",CIRCLEAREA);
                rt.addValue("BG_Area",BGCIRCLEAREA);
                rt.addValue("Intensity_"+channelNames[j],averages_r[i][j]);
                rt.addValue("Mean_Intensity_"+channelNames[j],averages_r[i][j] / CIRCLEAREA);
                rt.addValue("Intensity_"+channelNames[j]+"_bg",averages_r_dil[i][j]);
                rt.addValue("Mean_Intensity_"+channelNames[j]+"_bg",averages_r_dil[i][j] / BGCIRCLEAREA);
                rt.addValue("Relative_intensity_"+channelNames[j], (averages_r[i][j] / CIRCLEAREA) - (averages_r_dil[i][j]-averages_r[i][j]) / (BGCIRCLEAREA-CIRCLEAREA) );
            }
        }

        if(report){
            rt.show("Raw Results Sample "+sample);
        }


        ArrayList<ArrayList<Double>> hGrams = new ArrayList<ArrayList<Double>>(nCh-1);

        int ii = 0;

        for(int j=0;j<nCh;j++) {
            if(j!=refChNo) {

                hGrams.add(new ArrayList<Double>());

                for (int i = 0; i < cnt; i++) {
                    double value = rt.getValue("Relative_intensity_" + channelNames[j], i);
                    if (value < 0) {

                        hGrams.get(ii).add(value);
                        hGrams.get(ii).add(Math.abs(value));
                    }
                }
                ii++;
            }
        }

        ii=0;

        if(cutoffs==null) {

            cutoffs = new double[nCh];

            for (int i = 0; i < nCh; i++) {
                if (i != refChNo) {

                    Plot p = new Plot(channelNames[i], "x", "y");
                    p.addHistogram(hGramToArray(hGrams.get(ii)), getSizeCal.getFDbinsize(hGramToArray(hGrams.get(ii))) / 3.0);
                    float[] x_raw = p.getXValues();
                    float[] y_raw = p.getYValues();
                    double[] x = new double[x_raw.length];
                    double[] y = new double[y_raw.length];

                    double xMin = (double) x_raw[0];
                    double xMax = (double) x_raw[0];

                    for (int j = 0; j < x_raw.length; j++) {
                        x[j] = (double) x_raw[j];
                        y[j] = (double) y_raw[j];
                        if (x[j] > xMax) {
                            xMax = x[j];
                        }
                        if (x[j] < xMin) {
                            xMin = x[j];
                        }
                    }

                    CurveFitter cv = new CurveFitter(x, y);
                    cv.doFit(CurveFitter.GAUSSIAN);
                    double[] results = cv.getParams();

                    double[] xvalues = new double[100];
                    double[] yvalues = new double[100];

                    for (int j = 0; j < xvalues.length; j++) {

                        xvalues[j] = (double) (j * (xMax - xMin) / 100.0) + xMin;
                        yvalues[j] = results[0] + (results[1] - results[0]) * Math.exp(-Math.pow(xvalues[j] - results[2], 2) / (2 * Math.pow(results[3], 2)));

                    }
                    p.setColor(Color.RED);
                    p.add("line", xvalues, yvalues);

                    if (report) {
                        p.show();
                    }

                    cutoffs[i] = results[3] * 3;

                    ii++;
                }
            }

        }


        positives = new boolean[rt.size()][nCh];
        double[][] rawvalues = new double[nCh][rt.size()];
        ArrayList<Double> v = new ArrayList<Double>();
        medians = new double[nCh];

        for(int i=0;i<nCh;i++){
            v = new ArrayList<Double>();
            for(int j=0;j<rt.size();j++){

                double tval = rt.getValue("Relative_intensity_" + channelNames[i] ,j);
                if(tval > cutoffs[i]){
                positives[j][i] = true;
                v.add(tval);
                }

            }

            rawvalues[i] = hGramToArray(v);
            medians[i] = EVQ_sizing_analysis.getMedian(rawvalues[i]);

        }

        int NumberOfCat = 0;

        for(int i=0;i<nCh;i++){
            NumberOfCat += Math.pow(2,i);
        }

        int[] catCount = new int[NumberOfCat+1];
        int[] catCountch = new int[nCh];


        for(int i=0;i<positives.length;i++){

            int cat = 0;

            for(int j=0;j<positives[i].length;j++){
                int a = positives[i][ positives[i].length-j-1 ] ? 1 : 0;
                cat = cat +  (a * (int) Math.pow(2,j));

                if(positives[i][positives[i].length-j-1]){
                    catCountch[j]++;
                }


            }

            catCount[cat]++;


        }

        int catLength=0;

        for(int i=0;i<catCount.length;i++){
            String cat = Integer.toBinaryString(i);
            String catName = "";

            while(cat.length() < nCh ){
                cat = "0"+cat;
            }

            String[] catstr = cat.split("");

            if(catstr[refChNo].equals("1")){
                catLength++;
            }
        }

        thrs = cutoffs;
        ChNames = channelNames;
        CollumnNames = new String[(catLength + catCountch.length )*3];
        Values = new double[CollumnNames.length];
        int fillcounter =0;



        double conc = ((1000000000000.0 /(  (impCh[0].getWidth()*c.pixelWidth) * (impCh[0].getHeight()*c.pixelHeight) * dups * eots  )) * (totalVolume/sampleVolume) * dil) * positives.length;


        CollumnNames[fillcounter] = "Total";
        Values[fillcounter] = (double) positives.length;
        fillcounter++;

        CollumnNames[fillcounter] = "Total_pecr";
        Values[fillcounter] = (double) 100;
        fillcounter++;

        CollumnNames[fillcounter] = "Total_conc";
        Values[fillcounter] = conc;
        fillcounter++;


        for(int i=0;i<catCountch.length;i++){

            if(i != refChNo) {

                String catName = "+"+ channelNames[i];

                conc = ((1000000000000.0 / ((impCh[0].getWidth() * c.pixelWidth) * (impCh[0].getHeight() * c.pixelHeight) * eots * dups)) * (totalVolume / sampleVolume) * dil) * catCountch[catCountch.length - i - 1];


                CollumnNames[fillcounter] = catName;
                Values[fillcounter] = (double) catCountch[catCountch.length - i - 1];
                fillcounter++;

                CollumnNames[fillcounter] = catName+"_perc";
                Values[fillcounter] = (double) catCountch[catCountch.length - i - 1] / (double) positives.length * 100.0;
                fillcounter++;

                CollumnNames[fillcounter] = catName+"_conc";
                Values[fillcounter] = conc;
                fillcounter++;

            }
        }


        for(int i=0;i<catCount.length;i++){

            String cat = Integer.toBinaryString(i);
            String catName = "";

            while(cat.length() < nCh ){
                cat = "0"+cat;
            }

            String[] catstr = cat.split("");




            for(int j=0;j<catstr.length;j++){

                if(j!=refChNo && catstr[j].equals("1")){
                    catName = catName+"+"+channelNames[j];
                }
                if(j!=refChNo && catstr[j].equals("0")){
                    catName = catName+"-"+channelNames[j];

                }


            }



            conc = ((1000000000000.0 /(  (impCh[0].getWidth()*c.pixelWidth) * (impCh[0].getHeight()*c.pixelHeight) * eots  ))* (totalVolume/sampleVolume) * dil) * catCount[i];


            if(catstr[refChNo].equals("1") && !catName.equals("")) {

                CollumnNames[fillcounter] = catName;
                Values[fillcounter] = catCount[i];
                fillcounter++;

                CollumnNames[fillcounter] = catName+"_perc";
                Values[fillcounter] = (double)catCount[i]/(double)positives.length*100.0;
                fillcounter++;

                CollumnNames[fillcounter] = catName+"_conc";
                Values[fillcounter] = conc;
                fillcounter++;

            }

        }






        if(report) {
            for (int i = 0; i < nCh; i++) {
                if (i != refChNo) {
                    IJ.log("CutOff values used for " + channelNames[i] + " :" + cutoffs[i]);
                }
            }
        }


    }





    private double[] hGramToArray( ArrayList<Double> hGram){

        double[] out = new double[hGram.size()];

        for (int i = 0; i < hGram.size(); i++) {
            out[i] = hGram.get(i);
        }

        return out;
    }
}
