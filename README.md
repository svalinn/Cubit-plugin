Svalinn plugins and command extensions for Cubit
=================================================

The plugin has been tested and is confirmed to work with versions of Coreform
Cubit up to 2021.5 and versions of Cubit up to 15.8.

|              | Trelis 17.01       | Cubit 2021.4 | Cubit 2021.5 | Cubit 2021.11 |
|--------------|--------------------|--------------|--------------|--------------|
| Ubuntu 18.04 | [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_ubuntu-18.04_cubit_17.1.0.tgz)|              |              |              |
| Ubuntu 20.04 |                    |  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_ubuntu-20.04_cubit_2021.4.tgz)|  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_ubuntu-20.04_cubit_2021.5.tgz) | [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_ubuntu-20.04_cubit_2021.11.tgz) |
| Ubuntu 21.04 |                    |  [![version](https://img.shields.io/badge/version-0.2.2-orange)](https://github.com/svalinn/Cubit-plugin/releases/download/0.2.2/svalinn-plugin_ubuntu-21.04_cubit_2021.4.tgz)|  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_ubuntu-21.04_cubit_2021.5.tgz) | [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_ubuntu-21.04_cubit_2021.11.tgz) |
| Debian 10.10 |                    |              |  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_debian-10.10_cubit_2021.5.tgz) | [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_debian-10.10_cubit_2021.11.tgz) |
| Mac OS  |  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_macos_cubit_17.1.0.tgz) |  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_macos_cubit_2021.4.tgz) |  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_macos_cubit_2021.5.tgz) | [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn-plugin_macos_cubit_2021.11.tgz) |
|Windows  |  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn_plugin_windows_17.1.0.zip) |  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn_plugin_windows_2021.4.zip) |  [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn_plugin_windows_2021.5.zip) | [![version](https://img.shields.io/badge/version-0.2.3-green)](https://github.com/svalinn/Cubit-plugin/releases/download/v0.2.3/svalinn_plugin_windows_2021.11.zip) |

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

Installation with Coreform Cubit (as of version 2023.8)
=======================================================

The DAGMC plugin is now shipped with Coreform Cubit starting at version 2023.8, but it is not activated by default. Following the instructions corresponding to your operating system to enable the plugin.

On Windows and Linux, use the dropdown menu "Tools", then select "Plugins...". Add the bin/plugins directory via the diaglog box. Restart Cubit to activate the plugin.

On Mac, use the dropdown menu "Tools", then select plugins and select "Add". Next, in Finder, navigate to /Applications/Coreform-Cubit-2023.8.app/Contents/lib/plugins.
Drag the plugins folder from the finder and drop it in the Cubit window. Once you restart Cubit, the plugin will be activated. You can confirm this by typing "help dagmc" in the command line.


Manually Install The Plugin
===========================

The Plugin can be downloaded from the [Release](https://github.com/svalinn/Cubit-plugin/releases)
section of this repository. Operating system specific assets are created with
each release and can be download for various Operating systems (Ubuntu 18, 20
and MacOS 10.15).

After downloading the compressed tgz file the contents should be uncompress to
your Cubit folder.
<details>
    <summary><b>Expand</b> - Terminal commands for Ubuntu 20.04</summary>
    <pre><code class="language-html">
    wget https://github.com/svalinn/Cubit-plugin/releases/download/0.1.0/svalinn-plugin_ubuntu-20.04_cubit_2021.5.tgz
    sudo tar -xzvf svalinn-plugin_ubuntu-20.04_cubit_2021.5.tgz -C /opt/Coreform-Cubit-2021.5
    </code></pre>
</details>

Test the Plugin
===============

Load Cubit by running ```coreform_cubit``` from the command line. After the
Cubit GUI finishes loading, run the following command in the Cubit command line:

```
help dagmc
```

If the plugin was installed correctly you should see the following text appear.

```
Help for words: dagmc.

export dagmc <filename> [faceting_tolerance <faceting tolerance>] [length_tolerance <length tolerance>]
     [normal_tolerance <normal tolerance>] [make_watertight] [verbose]
     [fatal_on_curves]
```

If this output does not appear, then the plugin was not installed correctly.

Usage within Cubit
==================

The DAGMC plugin commands should now be accessable via the Cubit command line.

To view the available command line options for the MCNP importer, type
```help mcnp``` in the Cubit command line.

Similarly for the DAGMC exporter, type ```help dagmc``` in the Cubit command line.

Some sample files have been included in the [test_plugin](test_plugin)
directory of this repository for testing the command line usage.


Building the plugin
===================

Building the plugin from source is also [possible](README_dev.md) but only recommended for developers.


Updating the repository for new Cubit releases
==============================================

Add the Cubit new version number to the following three lists:

- [Linux yml file](https://github.com/svalinn/Cubit-plugin/blob/453a2903306a635dbaedb573521f96351a83ed6b/.github/workflows/unix_linux.yml#L35)

- [Mac yml file](https://github.com/svalinn/Cubit-plugin/blob/453a2903306a635dbaedb573521f96351a83ed6b/.github/workflows/unix_mac.yml#L37)

- [Windows file](https://github.com/svalinn/Cubit-plugin/blob/453a2903306a635dbaedb573521f96351a83ed6b/.github/workflows/windows.yml#L38)

Create a new ```elif``` entry in the ```Environment Variables``` stage of the yml files. Set the BASE variable equal to the full URL (with the hash) of the new Cubit release.

- [Linux yml file](https://github.com/svalinn/Cubit-plugin/blob/453a2903306a635dbaedb573521f96351a83ed6b/.github/workflows/unix_linux.yml#L75-L92)

- [Mac yml file](https://github.com/svalinn/Cubit-plugin/blob/453a2903306a635dbaedb573521f96351a83ed6b/.github/workflows/unix_mac.yml#L49-L66)

- [Windows file](https://github.com/svalinn/Cubit-plugin/blob/453a2903306a635dbaedb573521f96351a83ed6b/.github/workflows/windows.yml#L53-L67)

Add an entry for the new Cubit version in the ```scripts/unix_share_build.sh``` file

- [unix build sh](https://github.com/svalinn/Cubit-plugin/blob/453a2903306a635dbaedb573521f96351a83ed6b/scripts/unix_share_build.sh#L218)

Add new links in the table at the top of the README.md to contain a link to the new plugin files (tgz or zip files). These files won't exist at this stage but the links can be anticipated in advance (see existing links for the URL pattern).

Create a pull request to the develop branch with these changes to the .github/workflows CI yml files and README.md

Wait for the CI to run tests on the new code and the pull requested to be approved

Once the pull request is merged into develop then create a release from the develop branch, this will launch more CI that builds the plugins for various systems and versions. The links added to the README.md should now allow the newly build plugin to be downloaded.
