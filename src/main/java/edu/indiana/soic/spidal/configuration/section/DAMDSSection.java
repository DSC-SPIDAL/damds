package edu.indiana.soic.spidal.configuration.section;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

public class DAMDSSection {
    public DAMDSSection(String configurationFilePath) {
        Properties p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));
            distanceMatrixFile = p.getProperty("DistanceMatrixFile", "distance.bin");
            weightMatrixFile = p.getProperty("WeightMatrixFile", "distance.bin");
            labelFile = p.getProperty("LabelFile", "labels.txt");
            pointsFile = p.getProperty("PointsFile", "points.txt");
            timingFile = p.getProperty("TimingFile", "timings.txt");
            summaryFile = p.getProperty("SummaryFile","summary.txt");

            numberDataPoints = Integer.parseInt(p.getProperty("NumberDataPoints","-1"));
            targetDimension = Integer.parseInt(p.getProperty("TargetDimension","3"));
            threshold = Double.parseDouble(p.getProperty("Threshold", "0.000001"));
            alpha = Double.parseDouble(p.getProperty("Alpha", "0.95"));
            cgIter = Integer.parseInt(p.getProperty("CGIterations", "20"));
            cgErrorThreshold = Double.parseDouble(p.getProperty("CGErrorThreshold", "1"));
            isSammon = Boolean.parseBoolean(p.getProperty("IsSammon", "false"));

            isBigEndian = Boolean.parseBoolean(p.getProperty("IsBigEndian", "false"));
            isMemoryMapped = Boolean.parseBoolean(p.getProperty("IsMemoryMapped", "true"));
        } catch (IOException e) {
            throw new RuntimeException("IO exception occurred while reading configuration properties file", e);
        }
    }

    public String distanceMatrixFile;
    public String weightMatrixFile;
    public String labelFile;
    public String pointsFile;
    public String timingFile;
    public String summaryFile;

    public int numberDataPoints;
    public int targetDimension;
    public double threshold;
    public double alpha;
    public int cgIter;
    public double cgErrorThreshold;
    public boolean isSammon;

    public boolean isBigEndian;
    public boolean isMemoryMapped;
}


