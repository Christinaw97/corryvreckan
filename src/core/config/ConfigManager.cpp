/**
 * @file
 * @brief Implementation of config manager
 *
 * @copyright Copyright (c) 2017-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "ConfigManager.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "Configuration.hpp"
#include "core/config/exceptions.h"
#include "core/utils/log.h"

using namespace corryvreckan;

/**
 * @throws ConfigFileUnavailableError If the main configuration file cannot be accessed
 */
ConfigManager::ConfigManager(std::filesystem::path file_name,
                             std::initializer_list<std::string> global,
                             std::initializer_list<std::string> ignore) {
    // Check if the file exists
    std::ifstream file(file_name);
    if(!file || !std::filesystem::is_regular_file(file_name)) {
        throw ConfigFileUnavailableError(file_name);
    }

    // Convert main file to absolute path
    file_name = std::filesystem::canonical(file_name);
    LOG(TRACE) << "Reading main configuration";

    // Read the file
    ConfigReader reader(file, file_name);

    // Convert all global and ignored names to lower case and store them
    auto lowercase = [](const std::string& in) {
        std::string out(in);
        std::transform(out.begin(), out.end(), out.begin(), ::tolower);
        return out;
    };
    std::transform(global.begin(), global.end(), std::inserter(global_names_, global_names_.end()), lowercase);
    std::transform(ignore.begin(), ignore.end(), std::inserter(ignore_names_, ignore_names_.end()), lowercase);

    // Initialize global base configuration
    global_config_ = reader.getHeaderConfiguration();

    // Store all the configurations read
    for(auto& config : reader.getConfigurations()) {
        // Skip all ignored sections
        std::string config_name = config.getName();
        std::transform(config_name.begin(), config_name.end(), config_name.begin(), ::tolower);
        if(ignore_names_.find(config_name) != ignore_names_.end()) {
            continue;
        }

        // Merge all global section with the global config
        if(global_names_.find(config_name) != global_names_.end()) {
            global_config_.merge(config);
            continue;
        }

        module_configs_.push_back(config);
    }
}

void ConfigManager::parse_detectors() {
    // Reading detector file
    auto detector_file_names = global_config_.getPathArray("detectors_file", true);
    LOG(TRACE) << "Reading detector configurations";

    for(auto& detector_file_name : detector_file_names) {
        std::ifstream detector_file(detector_file_name);
        ConfigReader detector_reader(detector_file, detector_file_name);
        auto detector_configs = detector_reader.getConfigurations();
        detector_configs_.splice(detector_configs_.end(),
                                 std::list<Configuration>(detector_configs.begin(), detector_configs.end()));
    }
}

/**
 * The global configuration is the combination of all sections with a global header.
 */
Configuration& ConfigManager::getGlobalConfiguration() { return global_config_; }

/**
 * Load all extra options that should be added on top of the configuration in the file. The options loaded here are
 * automatically applied to the module instance when these are added later.
 */
bool ConfigManager::loadModuleOptions(const std::vector<std::string>& options) {
    bool optionsApplied = false;

    // Parse the options
    for(const auto& option : options) {
        module_option_parser_.parseOption(option);
    }

    // Apply global options
    optionsApplied = module_option_parser_.applyGlobalOptions(global_config_) || optionsApplied;

    // Apply module options
    for(auto& config : module_configs_) {
        optionsApplied = module_option_parser_.applyOptions(config.getName(), config) || optionsApplied;
    }

    return optionsApplied;
}

/**
 * Load all extra options that should be added on top of the detector configuration in the file. The options loaded here are
 * automatically applied to the detector instance when these are added later and will be taken into account when possibly
 * loading customized detector models.
 */
bool ConfigManager::loadDetectorOptions(const std::vector<std::string>& options) {
    bool optionsApplied = false;

    // Create the parser
    OptionParser detector_option_parser;

    // Parse the options
    for(const auto& option : options) {
        detector_option_parser.parseOption(option);
    }

    if(detector_configs_.empty()) {
        parse_detectors();
    }

    // Apply detector options
    for(auto& config : detector_configs_) {
        optionsApplied = detector_option_parser.applyOptions(config.getName(), config) || optionsApplied;
    }

    return optionsApplied;
}

/**
 * All special global and ignored sections are not included in the list of module configurations.
 */
std::list<Configuration>& ConfigManager::getModuleConfigurations() { return module_configs_; }
/**
 * The list of detector configurations is read from the configuration defined in 'detector_file'
 */
std::list<Configuration>& ConfigManager::getDetectorConfigurations() {
    if(detector_configs_.empty()) {
        parse_detectors();
    }
    return detector_configs_;
}

/**
 * @warning A previously stored configuration is directly invalidated if the same unique name is used again
 *
 * An instance configuration is a specialized configuration for a particular module instance. If a ModuleIdentifier already
 * exists an error is thrown.
 */
Configuration& ConfigManager::addInstanceConfiguration(const ModuleIdentifier& identifier, const Configuration& config) {
    // Check uniqueness
    if(instance_identifier_to_config_.find(identifier) != instance_identifier_to_config_.end()) {
        throw ModuleIdentifierAlreadyAddedError(identifier);
    }

    // Add configuration
    instance_configs_.push_back(config);
    Configuration& ret_config = instance_configs_.back();
    instance_identifier_to_config_[identifier] = --instance_configs_.end();

    // Add identifier key to config
    ret_config.set<std::string>("identifier", identifier.getIdentifier());

    // Apply instance options
    module_option_parser_.applyOptions(identifier.getUniqueName(), ret_config);
    return ret_config;
}

/**
 * The list of instance configurations can contain configurations with duplicate names, but the instance configuration is
 * guaranteed to have a configuration value 'identifier' that contains an unique identifier for every same config name
 */
std::list<Configuration>& ConfigManager::getInstanceConfigurations() { return instance_configs_; }

/**
 * An instance configuration might be dropped when not used (e.g. it is overwritten by another module instance afterwards)
 * We need to remove it from the instance configuration list to ensure dumping the config actually dumps only the instance
 * configurations that were used.
 */
void ConfigManager::dropInstanceConfiguration(const ModuleIdentifier& identifier) {
    // Remove config from instance configs and from instance identifier map
    if(instance_identifier_to_config_.find(identifier) != instance_identifier_to_config_.end()) {
        instance_configs_.erase(instance_identifier_to_config_[identifier]);
        instance_identifier_to_config_.erase(identifier);
    } else {
        throw ModuleIdentifierNotFoundError(identifier);
    }
}
