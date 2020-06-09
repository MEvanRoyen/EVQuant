package EVQ;


import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

/*
EVQuant FIJI Plugin
written by J.A. Slotman @ www.ErasmusOIC.nl
Based on FIJI macro's and excel templates By T.A. Hartjes

Version 1.0.2096 : Initial upload to github
*/

import java.awt.*;

public class getEOST {

    int[][] allPoints;
    int[][] centerPoints;
    double EOST;

    public getEOST(ImagePlus input, int threshold, boolean report){

        MyOwnFindmaxima mf = new MyOwnFindmaxima();
        GaussianBlur gb = new GaussianBlur();
        gb.showProgress(false);
        ImagePlus imp = input.duplicate();
        Calibration c = imp.getCalibration();

        ImageStack ims = imp.getImageStack();
        ImageProcessor[] ips = new ImageProcessor[ims.size()];
        Polygon[] polys = new Polygon[ims.size()];



        for(int i=1;i<=ims.size();i++){
            ips[i-1] = ims.getProcessor(i);
            gb.blurGaussian(ips[i - 1], 2, 2, 0.002);
            int num = (int) mf.findMaxima(ips[i-1], threshold, 5, true);
            polys[i-1] = mf.xyCoordinates;
            IJ.showStatus("Detecting Particles...");
            IJ.showProgress(i,ims.size());
        }

        allPoints = polyToCoordinates(polys);
        makeCategories();
        int[] centers = getCenters();

        ResultsTable rt = new ResultsTable();
        Analyzer.setResultsTable(rt);

        double allcount_mean = 0;
        double centercount_mean = 0;
        int slicenum = 0;


        for(int i=ims.size()/4;i<=ims.size()-(ims.size()/4);i++){

            int allcount = 0;
            int centercount = 0;

            for(int j=0;j<allPoints.length;j++){
                if(allPoints[j][2]==i){
                    allcount++;
                }
            }

            for(int j=0;j<centers.length;j++){
                if(centers[j] == i){
                    centercount++;
                }
            }

            if(report) {
                rt.incrementCounter();
                rt.addValue("Slice", i);
                rt.addValue("All_counts", allcount);
                rt.addValue("Center_counts", centercount);
            }

            allcount_mean += allcount;
            centercount_mean += centercount;
            slicenum++;

            IJ.showStatus("Analyzing Particles...");
            IJ.showProgress(i-ims.size()/4 , ims.size());

        }



        allcount_mean = allcount_mean / (double) slicenum;
        centercount_mean = centercount_mean / (double) slicenum;



        EOST = (allcount_mean/centercount_mean) * c.pixelDepth ;




        if(report) {
            rt.show("Results");

            ImageStack ims_result = new ImageStack(ims.getWidth(), ims.getHeight());
            for (int z = 1; z <= ims.size(); z++) {
                ImageProcessor ip_result = new ByteProcessor(ims.getWidth(), ims.getHeight());
                ims_result.addSlice(ip_result);
            }

            for (int i = 0; i < allPoints.length; i++) {
                ims_result.getProcessor(allPoints[i][2]).putPixel(allPoints[i][0]-1, allPoints[i][1]-1, allPoints[i][3]);
                ims_result.getProcessor(allPoints[i][2]).putPixel(allPoints[i][0]-1, allPoints[i][1], allPoints[i][3]);
                ims_result.getProcessor(allPoints[i][2]).putPixel(allPoints[i][0]-1, allPoints[i][1]+1, allPoints[i][3]);
                ims_result.getProcessor(allPoints[i][2]).putPixel(allPoints[i][0]+1, allPoints[i][1]-1, allPoints[i][3]);
                ims_result.getProcessor(allPoints[i][2]).putPixel(allPoints[i][0]+1, allPoints[i][1], allPoints[i][3]);
                ims_result.getProcessor(allPoints[i][2]).putPixel(allPoints[i][0]+1, allPoints[i][1]+1, allPoints[i][3]);
                ims_result.getProcessor(allPoints[i][2]).putPixel(allPoints[i][0], allPoints[i][1]-1, allPoints[i][3]);
                ims_result.getProcessor(allPoints[i][2]).putPixel(allPoints[i][0], allPoints[i][1], allPoints[i][3]);
                ims_result.getProcessor(allPoints[i][2]).putPixel(allPoints[i][0], allPoints[i][1]+1, allPoints[i][3]);

            }

            ImagePlus imp_result = new ImagePlus("Result", ims_result);
            imp_result.show();

            IJ.log("Calculated EOST: "+EOST);
        }

        IJ.showProgress(1);

    }

    public int[][] polyToCoordinates(Polygon[] polys){

        int n = 0;

        for(int i=0;i<polys.length;i++){

            n = n + polys[i].npoints;

        }

        int[][] allpoints = new int[n][4];
        int index = 0;


        for(int i=0;i<polys.length;i++){

            for(int j=0;j<polys[i].npoints;j++){
                allpoints[index][0] = polys[i].xpoints[j];
                allpoints[index][1] = polys[i].ypoints[j];
                allpoints[index][2] = i+1;
                allpoints[index][3] = -1;

                index++;
            }

        }

        return allpoints;


    }

    public void makeCategories(){

        int cat=1;

        for(int i=0;i<allPoints.length-1;i++){

            if(allPoints[i][3]==-1){
                allPoints[i][3]=cat;
                cat++;
            }

            for(int j=i+1;j<allPoints.length;j++){
                if( inDist(allPoints[i][0],allPoints[i][1],allPoints[j][0],allPoints[j][1],7) && Math.abs(allPoints[i][2]-allPoints[j][2])==1 ) {
                    allPoints[j][3] = Math.max(allPoints[i][3],allPoints[i][3]);
                }
            }
        }


    }


    public int[] getCenters(){

        int maxCat = 0;

        for(int[] i : allPoints){
            if(i[3]>maxCat){
                maxCat=i[3];
            }
        }

        int[] centerCat = new int[maxCat-1];
        centerPoints = new int[maxCat-1][4];

        for(int i=0;i<centerCat.length;i++){

            double x = 0;
            double y = 0;
            double z = 0;

            double count = 0;

            for(int[] j : allPoints){
                if(j[3]==i+1){
                    x += j[0];
                    y += j[1];
                    z += j[2];
                    count++;
                }
            }

            x = (x/count);
            y = (y/count);
            z = (z/count);

            centerCat[i] = (int) z;
            centerPoints[i][0] = (int) x;
            centerPoints[i][1] = (int) y;
            centerPoints[i][2] = centerCat[i];
            centerPoints[i][3] = i+1;
        }

        return centerCat;

    }


    public boolean inDist(int x1, int y1, int x2, int y2, int dist){
        double powDist = Math.pow(x1-x2,2)+ Math.pow(y1-y2,2);

        if(powDist <= Math.pow(dist,2)){
            return true;
        }else{
            return false;
        }
    }


}
