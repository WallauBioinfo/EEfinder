import pkg_resources

package_info = pkg_resources.get_distribution("eefinder")
__version__ = package_info.version
