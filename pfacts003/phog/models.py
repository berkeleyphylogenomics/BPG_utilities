# DON'T DELETE THIS FILE!

# In Django Parlance, an installed app (as listed in settings.INSTALLED_APPS),
# is not actually an installed app if it doesn't have models. Django expects
# each app to have its own models. If the models don't exist, the app isn't
# actually considered "installed", even if the app is listed in the
# INSTALLED_APPS. Therefore, either refactor the models in question to this
# file, or leave this file blank. Do not, however, delete this file.

# The consequences of not having an app "installed" may not be ovious at first
# (the website may still function without error). However, other things, like
# database synchronization and tests will not work properly. At the time of this
# writing, we discovered that the tests were not being run simply because this
# models.py file didn't exist.
