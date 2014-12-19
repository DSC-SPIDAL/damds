package edu.indiana.soic.spidal.configuration;

import edu.indiana.soic.spidal.configuration.section.DAMDSSection;

public class ConfigurationMgr {
    private String configurationFilePath;
    public DAMDSSection damdsSection;

    public ConfigurationMgr(String configurationFilePath) {
        this.configurationFilePath = configurationFilePath;
        damdsSection = new DAMDSSection(configurationFilePath);


    }

    public static ConfigurationMgr LoadConfiguration(String configurationFilePath){
        // TODO - Fix configuration management
        return new ConfigurationMgr(configurationFilePath);
    }
}
