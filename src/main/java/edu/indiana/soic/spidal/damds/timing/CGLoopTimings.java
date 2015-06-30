package edu.indiana.soic.spidal.damds.timing;

import com.google.common.base.Stopwatch;

import java.util.concurrent.TimeUnit;

public class CGLoopTimings {
    public static enum TimingTask{
        MM, INNER_PROD_PAP, INNER_PROD_R
    }

    private static Stopwatch timerMM = Stopwatch.createUnstarted();
    private static Stopwatch timerInnerProdPAP = Stopwatch.createUnstarted();
    private static Stopwatch timerInnerProdR = Stopwatch.createUnstarted();

    private static long tMM;
    private static long tInnerProdPAP;
    private static long tInnerProdR;

    private static long countMM;
    private static long countInnerProdPAP;
    private static long countInnerProdR;

    public static void startTiming(TimingTask task){
        switch (task){
            case MM:
                timerMM.start();
                ++countMM;
                break;
            case INNER_PROD_PAP:
                timerInnerProdPAP.start();
                ++countInnerProdPAP;
                break;
            case INNER_PROD_R:
                timerInnerProdR.start();
                ++countInnerProdR;
                break;
        }
    }

    public static void endTiming(TimingTask task){
        switch (task){
            case MM:
                timerMM.stop();
                tMM += timerMM.elapsed(TimeUnit.MILLISECONDS);
                timerMM.reset();
                break;
            case INNER_PROD_PAP:
                timerInnerProdPAP.stop();
                tInnerProdPAP += timerInnerProdPAP.elapsed(TimeUnit.MILLISECONDS);
                timerInnerProdPAP.reset();
                break;
            case INNER_PROD_R:
                timerInnerProdR.stop();
                tInnerProdR += timerInnerProdR.elapsed(TimeUnit.MILLISECONDS);
                timerInnerProdR.reset();
                break;
        }
    }

    public static double getTotalTime(TimingTask task){
        switch (task){
            case MM:
                return tMM;
            case INNER_PROD_PAP:
                return tInnerProdPAP;
            case INNER_PROD_R:
                return tInnerProdR;
        }
        return  0.0;
    }

    public static double getAverageTime(TimingTask task){
        switch (task){
            case MM:
                return tMM *1.0/ countMM;
            case INNER_PROD_PAP:
                return tInnerProdPAP *1.0/ countInnerProdPAP;
            case INNER_PROD_R:
                return tInnerProdR *1.0/ countInnerProdR;
        }
        return  0.0;
    }

}
