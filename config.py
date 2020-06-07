import os


class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'
    ALLOWED_IMAGE_EXTENSIONS = [".VCF"]
    DATA_FOLDER = "./custom_panel"
    DEPENDENCIES_FOLDER = "./dependencies"
