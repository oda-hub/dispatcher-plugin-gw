Gravitational Waves MMODA plugin
========

What's the license?
-------------------

GW plugin is distributed under the terms of The MIT License.

Who's responsible?
-------------------
Denys Savchenko

Configuration for deployment
----------------------------
- copy the `conf_file` from `dispatcher_plugin_gw/config_dir/data_server_conf.yml' and place in given directory
- set the environment variable `CDCI_GW_PLUGIN_CONF_FILE` to the path of the file `conf_file` 
- edit the in `conf_file` the key:
    - `data_server_url:`  