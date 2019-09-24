import os


if __name__ == "__main__":
    aquapi_username = os.environ.get('AQUAPI_USERNAME', None)
    aquapi_password = os.environ.get('AQUAPI_PASSWORD', None)
    aquapi_url = os.environ.get('AQUAPI_URL', None)
    release_python = os.environ.get('RELEASE_PYTHON', 'False')

    if release_python == 'True':
        os.system(f'devpi use {aquapi_url}')
        os.system(f'devpi login {aquapi_username} --password {aquapi_password}')
        os.system('devpi upload --no-vcs --formats sdist')
