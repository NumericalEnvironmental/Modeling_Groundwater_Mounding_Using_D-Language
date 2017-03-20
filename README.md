# Modeling_Groundwater_Mounding_Using_D-Language

![Preview](https://numericalenvironmental.files.wordpress.com/2017/03/model-comparison.jpg?w=816)

This is a relatively short program, written and compiled in D-language (https://dlang.org/), which simulates time-dependent mounding of groundwater in an unconfined aquifer underneath recharge basin(s), subject to pumping from nearby well(s). The mounding solution is based on the Hantush (1967) model; the well model is based on the familiar Theis (1935) solution. The models are solved separately (with special numerical routines needed for both) and then solutions are superimposed, along with a regional groundwater gradient if warranted.

The source code is included here for perusal, along with an executable that has been compiled to run under Windows. You will need to compile the model yourself if you are using a different platform. Required text input files include the following:

* aquifer.txt –basic aquifer properties (hydraulic conductivity, specific yield, starting saturated thickness, gradient)
* model.txt – model parameters (gridding, time stepping, whether or not to use particle tracking, etc.)
* basins.txt – recharge basin information, one basin per row (centroid location, length and width, rotation angle with respect to due east)
* wells.txt – well names and locations
* files for individual recharge basins and for wells, as names in the basins.txt and wells.txt files, with flux history tables (e.g., recharge rates or pumping rates as a function of time, using step functions)
* monitor_locs.txt – name sand locations of monitoring points (i.e., locations where continuous water level time series are recorded)
* particles_init.txt – initial positions of particles (for particle tracking, if this option is selected in the model.txt file)

More discussion is provided, along with an example application problem, here: https://numericalenvironmental.wordpress.com/2017/03/20/coding-an-analytical-solution-for-groundwater-banking-with-superposition/. The example problem includes a comparison to an equivalent simulation using MODFLOW; its input files (implemented using the USGS’s ModelMuse MODFLOW pre/post-processor) are included here as well under the MODFLOW subfolder.

Email me with questions at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

