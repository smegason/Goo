Developer's Guide
=================

We welcome everyone who is interested in creating 3D models of biological stem cells on Blender!

The developer's guide is to help contributors on how to start working and contributing 
on Goo. If you are interested in join, please contact Sean Megason on his email address. sean.megason(AT)hms.harvard.edu.

Git/Github
----------

The first thing you need to do is to fork Goo repository and then clone the forked repository into your local machine. This 
Will give you the freedom to test and make changes to the project without affecting the original repository. 

Make sure you have git installed on your local machine.  If you do not have it install on your system, refer to this guide: https://github.com/git-guides/install-git

Fork & clone
~~~~~~~~~~~~~~

If you are new to Git/Github, here is detailed info on how to fork and clone a repository from Github:
https://docs.github.com/en/get-started/quickstart/fork-a-repo

Push
~~~~~~

Once you make final changes and the scripts run with no issues, you can push them into the main project. Here_ is an instruction on how to push using git.

.. _Here: https://docs.github.com/en/get-started/importing-your-projects-to-github/importing-source-code-to-github/adding-locally-hosted-code-to-github 

Blender 
-------

Blender is an open source 3D computer software. It is free to download_. Blender has a huge community and contributors all over the world. You may find Python API_ for Blender very useful.  

.. _download: https://www.blender.org/download/
.. _API: https://docs.blender.org/api/current/index.html

Once you have downloaded and installed Blender, launch Blender and select *General* for the new file. **Scripting** tab allows you to interact with Blender using Python language. By clicking **scripting**, a list of windows will open:

- *3D Viewport* allows you to interact with 3D objects.
- *Python Interactive Console* allows you to interact with Blender using Python language. You can think of it as a command window/terminal of Blender. 
- *Info Console Menu* or *Report console* creates a log in Python of the executed actions on Blender. In another words, it translates the clicks you made in Blender into code. This is very useful especially when you try to automate procedures of actions on *scripting*.    

.. image:: ../img/basic_blender.png
  :width: 800
    

    