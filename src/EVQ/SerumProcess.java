package EVQ;


import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.ChannelSplitter;
import ij.plugin.RGBStackMerge;
import ij.process.ByteProcessor;

/*
EVQuant FIJI Plugin
written by J.A. Slotman @ www.ErasmusOIC.nl
Based on FIJI macro's and excel templates By T.A. Hartjes
*/

public class SerumProcess {

    ImagePlus imp;

    public SerumProcess(ImagePlus imp_){

        ImagePlus[] impsRaw = ChannelSplitter.split(imp_);
        String[] ChNames = new String[impsRaw.length];
        String[] label = {"Max", "Average"};
        int method;
        boolean[] SelectedCh = new boolean[impsRaw.length];

        for(int i=0; i<impsRaw.length; i++){
            ChNames[i] = "Channel "+(i+1);
        }

        GenericDialog gd = new GenericDialog("EVQuant Serum and Plasma analysis preprocessing");
        gd.addMessage("Select channels used for counting Vesicles don't include Rhodamine Channel");
        gd.addCheckboxGroup(1,impsRaw.length,ChNames,SelectedCh);
        gd.addChoice("Method",label,"Max");
        gd.showDialog();

        if(gd.wasCanceled()){
            return;
        }
        int numCh = 0;
        for(int i=0;i<SelectedCh.length;i++) {
            SelectedCh[i] = gd.getNextBoolean();
            if(SelectedCh[i]){
                numCh++;
            }
        }

        method = gd.getNextChoiceIndex();

        ImagePlus[] imps = new ImagePlus[numCh];


        int cnt = 0;
        for(int i=0;i<SelectedCh.length;i++){

            if(SelectedCh[i]){
                imps[cnt] = impsRaw[i];
                cnt++;
            }

        }

        ImagePlus MergedImp = new ImagePlus("Merged");
        ImageStack MergedIms = new ImageStack(imps[0].getWidth(),imps[0].getHeight(),imps[0].getImageStackSize());
        int[] MergedPixels = new int[imps[0].getWidth()*imps[0].getHeight()];
        byte[] MergedPixelsB = new byte[imps[0].getWidth()*imps[0].getHeight()];
        byte[] OriPixels = (byte[]) imps[0].getImageStack().getProcessor(1).getPixels();

        for(int j=0;j<imps[0].getImageStack().size();j++){

            for(int k=0;k<OriPixels.length;k++) {
                MergedPixels[k] = 0;
            }

            for(int i=0;i<imps.length;i++){
                OriPixels = (byte[]) imps[i].getImageStack().getProcessor(j+1).getPixels();
                for(int k=0;k<OriPixels.length;k++){
                    if(method==0) {
                        if (toUnsigned(OriPixels[k]) > MergedPixels[k]) {
                            MergedPixels[k] = toUnsigned(OriPixels[k]);
                        }
                    }
                    if(method==1) {
                        MergedPixels[k] += toUnsigned(OriPixels[k]);
                    }
                }


            }
            if(method==1){
                for(int k=0;k<OriPixels.length;k++) {
                    MergedPixels[k] /= imps.length;
                }
            }

            for(int k=0;k<OriPixels.length;k++) {
                MergedPixelsB[k] = toSigned(MergedPixels[k]);
            }

            ByteProcessor ip = new ByteProcessor(imps[0].getWidth(),imps[0].getHeight());
            ip.setPixels(MergedPixelsB.clone());
            MergedIms.setProcessor(ip,j+1);

        }

        MergedImp.setStack(MergedIms);


        ImagePlus[] OutputImp = new ImagePlus[imps.length+1];

        for(int i=0;i<OutputImp.length;i++){
            if(i==0){
                OutputImp[i] = MergedImp;
            }else{
                OutputImp[i] = imps[i-1];
            }
        }


        imp = RGBStackMerge.mergeChannels(OutputImp,false);









    }

    public ImagePlus getImp(){
        return imp;
    }

    public int toUnsigned(byte b){
        int out = b < 0 ? b+256 : b;
        return out;
    }

    public byte toSigned(int i){
        byte out = i < 128 ? (byte) i : (byte)(i-256);
        return out;
    }
}
