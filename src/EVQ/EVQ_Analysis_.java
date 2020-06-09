package EVQ;

import ij.*;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.PlugIn;
import ij.plugin.filter.GaussianBlur;
import ij.process.ImageProcessor;

/*
EVQuant FIJI Plugin
written by J.A. Slotman @ www.ErasmusOIC.nl
Based on FIJI macro's and excel templates By T.A. Hartjes
*/

public class EVQ_Analysis_ implements PlugIn {



    public void run(String args){

        boolean immuno = args.equals("immuno");
        boolean serum = args.equals("serum");


        double eots = Prefs.get("EVQ.EOST",3.0);
        int threshold = (int) Prefs.get("EVQ.Threshold",4);
        int dup = (int) Prefs.get("EVQ.dup",10);
        boolean edge = false;
        int edgesize = (int) Prefs.get("EVQ.bordersize",0);
        boolean raw = false;
        int outputType = 5;
        int samplevolume = (int) Prefs.get("EVQ.samplevolume",330);
        String[] channelNames = new String[1];

        ImagePlus imp_orig = WindowManager.getCurrentImage();

        if(serum){

            immuno = true;
            SerumProcess sp = new SerumProcess(imp_orig);
            imp_orig = sp.getImp();

        }

        if(imp_orig==null){
            IJ.showMessage("Error","No image is open");
            return;
        }

        ImagePlus imp = imp_orig.duplicate();
        ImagePlus[] impCh = new ImagePlus[1];

        if(immuno){
            impCh = ChannelSplitter.split(imp);
            imp = impCh[0];

            channelNames = new String[impCh.length];
        }


        ImageStack ims = imp.getStack();
        Calibration c = imp.getCalibration();


        GenericDialog gd = new GenericDialog("EVQuant Analysis setup");
        if(ims.getSize()!=1) {
            gd.addMessage("Loaded dataset " + imp.getTitle() + " with " + ims.getSize() + " samples");
        }else{
            gd.addMessage("Loaded dataset " + imp.getTitle() + " with " + ims.getSize() + " sample");
        }
        gd.addNumericField("Threshold",threshold,0);
        gd.addNumericField("EOST",eots,3);
        gd.addNumericField("Number of images per sample",dup,0);
        gd.addNumericField("Total volume in each well",samplevolume,0);
        if(immuno){
            for(int i=0;i<impCh.length;i++){
                channelNames[i] = "Ch"+(i+1);
                gd.addStringField("Channel "+(i+1)+" name",channelNames[i]);
            }
            gd.addRadioButtonGroup("Select Generic Dye channel",channelNames,1,impCh.length,channelNames[0]);
            gd.addCheckbox("Exclude on border",edge);
            gd.addToSameRow();
            gd.addNumericField("Size of border",edgesize,0);
        }
        gd.addCheckbox("Show Raw Results",raw);

        gd.showDialog();

        if(gd.wasCanceled()){
            return;
        }



        threshold = (int) gd.getNextNumber();
        eots = gd.getNextNumber();
        dup = (int) gd.getNextNumber();
        samplevolume = (int) gd.getNextNumber();
        String refCh;
        int refChNo = 0;
        if(immuno){
            for(int i=0;i<impCh.length;i++){
                channelNames[i] = gd.getNextString();
            }
            refCh = gd.getNextRadioButton();
            for(int i=0;i<impCh.length;i++){
                if(refCh.equals("Ch"+(i+1))){
                    refChNo = i;
                }
            }
            edge = gd.getNextBoolean();
            edgesize = (int) gd.getNextNumber();
        }



        raw = gd.getNextBoolean();

        if(!edge){
            edgesize=0;
        }

        if(edgesize < 0){
            IJ.showMessage("Error","Border cannot be negative, a border of 0 will be used");
            edgesize = 0;
        }


        if(immuno) {
            imp = impCh[refChNo];
        }

        Prefs.set("EVQ.Threshold",threshold);
        Prefs.set("EVQ.EOST",eots);
        Prefs.set("EVQ.dup",dup);
        Prefs.set("EVQ.samplevolume",samplevolume);
        Prefs.set("EVQ.bordersize",edgesize);



        int samplenum = ims.size()/dup;
        if(ims.size()%dup != 0){
            samplenum = samplenum + 1;
        }

        int[] sampledups = new int[samplenum];
        double[] sampledils = new double[samplenum];
        double[] samplevols = new double[samplenum];
        String[] samplenames = new String[samplenum];
        boolean[] blancos = new boolean[samplenum];
        int sampletot = 0;


        for(int i=0;i<samplenames.length;i++){
            samplenames[i] = "Sample_"+(i+1);
        }



        int j=0;
        int k=0;
        while(j<samplenames.length) {

            GenericDialog gd2 = new GenericDialog("EVQuant Analysis Sample setup");

            for (int i = 0; i < 10 && j<samplenames.length; i++) {
                gd2.addStringField(samplenames[j], samplenames[j]);
                gd2.addToSameRow();
                gd2.addNumericField("number of images", dup, 0);
                gd2.addToSameRow();
                gd2.addNumericField("sample volume (ul)",40,2);
                gd2.addToSameRow();
                gd2.addNumericField("dilution", 1, 3);
                gd2.addToSameRow();
                gd2.addCheckbox("Negative Control (NC)",false);
                j++;
            }

            gd2.showDialog();

            if (gd2.wasCanceled()) {
                return;
            }

            for (int i = 0; i < 10 && k<samplenames.length; i++) {

                samplenames[k] = gd2.getNextString();
                sampledups[k] = (int) gd2.getNextNumber();
                samplevols[k] = gd2.getNextNumber();
                sampledils[k] = gd2.getNextNumber();
                blancos[k] = gd2.getNextBoolean();
                sampletot += sampledups[k];
                k++;
            }


        }

        if(sampletot > ims.size()){
            IJ.showMessage("Error","More samples defined than slices present");
            return;
        }
        if(sampletot < ims.size()){
            IJ.log("(Warning) Not all slices are analyzed, because there were more slices than defined for the samples");
        }


        //analysis part

        MyOwnFindmaxima mf = new MyOwnFindmaxima();
        GaussianBlur gb = new GaussianBlur();
        gb.showProgress(false);

        int[] allvalues = new int[ims.size()];
        int[] samplevalues = new int[samplenum];
        double[] samplevaluesbg = new double[samplenum];

        int sample_index = 0;
        int sample_num_index = 0;
        if(!immuno) {
            for (int i = 1; i <= sampletot; i++) {

                ImageProcessor ip = ims.getProcessor(i);
                gb.blurGaussian(ip, 2, 2, 0.002);
                allvalues[i - 1] = (int) mf.findMaxima(ip, threshold, outputType, true);

                samplevalues[sample_index] = samplevalues[sample_index] + allvalues[i - 1];
                sample_num_index++;

                if (sample_num_index >= sampledups[sample_index]) {
                    sample_num_index = 0;
                    sample_index++;
                }

                IJ.showProgress(i, sampletot - 1);

            }
        }

        EVQ_imunno[] im = new EVQ_imunno[samplenum];

        if(immuno){
            ImageStack[] imsCh = new ImageStack[impCh.length];
            ImageStack[] imsCh_empty = new ImageStack[impCh.length];
            ImagePlus[] impCh_temp = new ImagePlus[impCh.length];

            for(int i=0;i<impCh.length;i++){
                imsCh_empty[i] = new ImageStack(impCh[i].getWidth(),impCh[i].getHeight());

            }

            imsCh = imsCh_empty;

            for (int i = 1; i <= sampletot; i++) {

                for(int ii=0;ii<impCh.length;ii++){
                    imsCh[ii].addSlice(impCh[ii].getImageStack().getProcessor(i));
                }
                sample_num_index++;


                if(sample_num_index >= sampledups[sample_index]){
                    for(int ii=0;ii<impCh.length;ii++){
                        impCh_temp[ii] = new ImagePlus("",imsCh[ii]);

                    }

                    im[sample_index] = new EVQ_imunno();
                    im[sample_index].immuno(impCh_temp, c,  refChNo,  threshold, raw ,  channelNames,eots, samplevolume,sampledups[sample_index],  samplevols[sample_index],  sampledils[sample_index], null,samplenames[sample_index],edgesize );
                    samplevalues[sample_index] = (int) im[sample_index].Values[0];

                    for(int iii=0;iii<impCh.length;iii++){
                        imsCh[iii] = new ImageStack(impCh[iii].getWidth(),impCh[iii].getHeight());
                    }

                    impCh_temp = new ImagePlus[impCh.length];

                    sample_index++;
                    sample_num_index = 0;
                }
            }

        }

        double blancoValue = 0;
        double blancoSlices = 0;
        int blancocnt = 0;

        double[] blancoValues = null;
        if(immuno){
            blancoValues = new double[im[0].Values.length];
        }


        for(int i=0;i<samplenames.length;i++){

            if(blancos[i]){
                blancoValue += samplevalues[i];
                blancoSlices += sampledups[i];
                blancocnt++;

                if(immuno){
                    for(int ii=0;ii<im[i].Values.length;ii++){
                        blancoValues[ii] += im[i].Values[ii];
                    }
                }
            }


        }
        if(blancocnt!=0) {
            blancoValue = blancoValue / blancoSlices;

            if(immuno){
                for(int ii=0;ii<im[0].Values.length;ii++){
                    blancoValues[ii] /= blancoSlices;
                }
            }
        }

        ResultsTable rt = new ResultsTable();

        for(int i=0;i<samplenames.length;i++){

            String blanco = "";

            samplevaluesbg[i] = samplevalues[i] - (blancoValue*sampledups[i]);

            if(blancocnt!=0 && !blancos[i]){

                blanco = "";
            }else{

                blanco = "NC";
            }




            rt.incrementCounter();
            rt.addValue("Sample name",samplenames[i]);
            rt.addValue("Sample volume",samplevols[i]);
            rt.addValue("Sample dilution",sampledils[i]);
            rt.addValue("Number of particles", samplevalues[i]);



            double conc = ((1000000000000.0 /(  ((ims.getWidth()-(edgesize*2))*c.pixelWidth) * ((ims.getHeight()-(2*edgesize))*c.pixelHeight) * eots * sampledups[i] )) * ((double)samplevolume/samplevols[i]) * sampledils[i]) * samplevaluesbg[i];

            if (blancocnt != 0) {
                rt.addValue("Number of particles (corrected for NC)", samplevaluesbg[i]);
                rt.addValue("Concentration particles per ml (corrected for NC)", conc);
                rt.addValue("Negative control (NC)", blanco);
            } else {
                rt.addValue("Concentration particles per ml", conc);
            }

            if(immuno){
                for(int ii=0;ii<im[i].Values.length;ii++){
                    if(ii%3==0) {
                        im[i].Values[ii] = im[i].Values[ii]-(blancoValues[ii]*sampledups[i]);
                        rt.addValue(im[i].CollumnNames[ii],im[i].Values[ii]);
                        if(im[i].Values[0]>=1){
                            im[i].Values[ii + 1] = im[i].Values[ii] / im[i].Values[0] * 100;
                        }else{
                            im[i].Values[ii + 1] = -1;
                        }
                        conc = ((1000000000000.0 /(  ((ims.getWidth()-(edgesize*2))*c.pixelWidth) * ((ims.getHeight()-(2*edgesize))*c.pixelHeight) * eots * sampledups[i] )) * ((double)samplevolume/samplevols[i]) * sampledils[i]) * im[i].Values[ii];
                        im[i].Values[ii+2] = conc;
                    }
                }

                for(int ii=0;ii<im[i].Values.length;ii++){
                    if(ii%3==1) {
                        if(im[i].Values[ii]==-1){
                            rt.addValue(im[i].CollumnNames[ii], "No particles");
                        }else {
                            rt.addValue(im[i].CollumnNames[ii], im[i].Values[ii]);
                        }
                    }
                }
                for(int ii=0;ii<im[i].Values.length;ii++){
                    if(ii%3==2) {
                        rt.addValue(im[i].CollumnNames[ii], im[i].Values[ii]);
                    }
                }

                for(int ii=0;ii<im[i].thrs.length;ii++){
                    if(ii!=refChNo) {
                        rt.addValue("Threshold " + im[i].ChNames[ii], im[i].thrs[ii]);
                    }
                }

                for(int ii=0;ii<im[i].medians.length;ii++){
                    if(ii!=refChNo) {
                        rt.addValue("Median " + im[i].ChNames[ii], im[i].medians[ii]);
                    }
                }




            }




        }

        rt.show("Results");

        if(raw && !immuno) {

            ResultsTable rt_raw = new ResultsTable();

            for (int i = 0; i < ims.size(); i++) {
                rt_raw.incrementCounter();
                rt_raw.addValue("Image Number", i);
                rt_raw.addValue("Number of Vesicles", allvalues[i]);
            }

            rt_raw.show("Raw Results");

        }



    }
}
