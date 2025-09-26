# Contributing to the documentation 

To documentation is generated with sphinx. 
The source code to these webpages is a set of markdown documents that can be found under
`gw_analysis_tools/docs_sphinx`
(on the github repo, they are [here](https://github.com/scottperkins/gw_analysis_tools/tree/master/docs_sphinx)). 
If you ever have an issue you had to hack your way through (e.g. getting the code to work on windows), 
please consider adding your solution to the documentation! 

# Contributing to the project

Contributions to this source code would be very much appreciated (both corrections or additions)!
We use branches and pull requests to facilitate this kind of work, so please follow the follow steps for contributing to this project:

0. Optional - For corrections, file a bug report first so people know someone's working on it and are aware of the issue

1. Pull the latest code from the repository (probably from the master branch):
```bash
git pull
```

2. Create a new branch and make it available on the online repository with the following command (where \<feature-branch\> can be anything you'd like to call it)
```bash
git checkout -b <feature-branch>
git push origin <feature-branch>
```

3. Next, make your corrections to the code

4. Stage the changes 
```bash
git add <modified/added files>
git commit -m "A very clear and helpful message"
```

5. Push the changes to the online repository
```bash
git push origin <feature-branch>
```

6. Make a pull request through the online GUI. This will start a discuss and verification process before the code is eventually merged by a moderator.
