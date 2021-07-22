Svalinn plugins and command extensions for Cubit
=================================================

The plugin has been tested and is confirmed to work with versions of Corefoam
Cubit up to 2021.5 and versions of Cubit up to 15.8.

Install Cubit
==============

Cubit or Corefoam Cubit can be installed by following the instructions on the
relevant website. Cubit is available for Sandia and US National labs users
while Coreform Cubit is available via comerical licenses or freeley available
for students and hobbyists as Cubit Learn.

- [Download Cubit](https://cubit.sandia.gov/downloads.html) from Sandia
- [Download Corefoam Cubit](https://coreform.com/products/downloads/) from Corefoam

<details>
    <summary><b>Expand</b> - Ubuntu terminal commands for Corefoam Cubit</summary>
    <pre><code class="language-html">
    sudo apt update
    sudo apt-get install wget
    wget -O coreform-cubit-2021.5.deb https://f002.backblazeb2.com/file/cubit-downloads/Coreform-Cubit/Releases/Linux/Coreform-Cubit-2021.5%2B15962_5043ef39-Lin64.deb
    sudo dpkg -i coreform-cubit-2021.5.deb 
    </code></pre>
</details>


Install Plugin
==============

The Plugin can be downloaded from the [Release](https://github.com/svalinn/Cubit-plugin/releases)
section of this repository. Operating system specific assets are created with
each release and can be download for various Operating systems (Ubuntu 18, 20
and MacOS 10.15).

After downloading the compressed tgz file the contents should be uncompress to
your Cubit folder.
<details>
    <summary><b>Expand</b> - Ubuntu terminal commands for Ubuntu 20.04</summary>
    <pre><code class="language-html">
    wget https://github.com/svalinn/Cubit-plugin/releases/download/0.1.0/svalinn-plugin_ubuntu-20.04_cubit_2021.5.tgz
    sudo tar -xzvf svalinn-plugin_ubuntu-20.04_cubit_2021.5.tgz -C /opt/Coreform-Cubit-2021.5
    </code></pre>
</details>


Test the Plugin
===============

Run ```coreform_cubit``` from the command line. If the plugin was installed
correctly, after the Cubit GUI finishes loading, the following output should
appear in the Cubit command line:

```
Loaded Svalinn plugin.
-- DAGMC export command available.
-- iGeom_test command available.
-- MCNP import command available.
Journaled Command: undo on

Cubit>
```

If this output does not appear, then the plugin was not installed correctly.


Usage within Cubit
==================

The DAGMC plugin commands should now be accessable via the Cubit command line.

To view the available command line options for the MCNP importer, type
```help mcnp``` in the Cubit command line.

Similarly for the DAGMC exporter, type ```help dagmc``` in the Cubit command line.

Some sample files have been included in the [test_plugin](https://github.com/svalinn/Cubit-plugin/tree/develop/test_plugin)
directory of this repository for testing the command line usage.


Building the plugin
===================

Building the plugin from source is also [possible](https://github.com/svalinn/Cubit-plugin/tree/develop/README_dev.md) but only recommended for developers.
