##
##  dshadow spec -- RPM Package Specification
##  This is a base SPEC file for dshadow. This file may be modified 
##  and written to a new file, if the rpms are built from the build_rom.sh
##  script in order to modify the Version field.
##
##  Permission to use, copy, modify, and distribute this software for
##  any purpose with or without fee is hereby granted, provided that
##  the above copyright notice and this permission notice appear in all
##  copies.
##
##  THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
##  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
##  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
##  IN NO EVENT SHALL THE AUTHORS AND COPYRIGHT HOLDERS AND THEIR
##  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
##  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
##  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
##  USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
##  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
##  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
##  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
##  SUCH DAMAGE.
##

#   package information
Name:         dshadow
Summary:      DShadow - A specialized fork of DaVinci
URL:          http://www.beneaththewaves.net/
Vendor:       Beneath the Waves
Packager:     Ben Lincoln <Use_the_Contact_form_on_the_website@beneaththewaves.net>  
Distribution: CentOS 5 (MSFF)
Group:        Applications/Science
License:      GPLv2
Version:      2.09
Release:      1

#   list of sources
Source:      http://www.beneaththewaves.net/Downloads/%{name}-%{version}.tar.bz2

#   build information
BuildPreReq:  automake, autoconf, make, gcc, gcc-c++, hdf5, hdf5-devel, cfitsio, cfitsio-devel, readline, readline-devel, zlib, zlib-devel, curl-devel	
BuildRoot: %{_tmppath}/%{name}-%{version}-buildroot
Requires:	gnuplot, hdf5, hdf5-devel, cfitsio, readline, zlib
#Provides:     
#Obsoletes:    
#Conflicts:    
Prefix: %_prefix

%description
	DaVinci's Shadow is a fork of ASU's DaVinci with some additional features.

%prep

%setup -q

%build
#Iomedley
cd iomedley
./configure --prefix=%{_prefix} --libdir=%{_libdir}  
make 
cd ../

#Davinci
./configure --prefix=%{_prefix} --libdir=%{_libdir} --disable-libisis   --with-viewer=/usr/bin/display \
--with-modpath=%{_libdir}/%{name} --with-help=%{_datadir}/%{name}/docs/dv.gih


make


%define __spec_install_post %{nil}
%define debug_package %{nil}

%install
rm -rf $RPM_BUILD_ROOT
mkdir -p -m 755 $RPM_BUILD_ROOT%{_bindir}
mkdir -p -m 755 $RPM_BUILD_ROOT%{_libdir}
make install DESTDIR=$RPM_BUILD_ROOT

%clean
rm -rf $RPM_BUILD_ROOT

%files
 %defattr(-,root,root)
 %{_bindir}/%{name}
 %{_libdir}/libdshadow*
 %{_libdir}/libiomedley*
 %{_libdir}/%{name}*
 %{_datadir}/%{name}*
 %{_includedir}/%{name}*





%post
##To avoid SELINUX  security message (maybe there is a better solution)
chcon -f -t textrel_shlib_t %{_libdir}/libdshadow* > /dev/null 2>&1 || /bin/true


 

%changelog
* Tue Jul 10 2007 Betim Deva <betim@asu.edu> 1.6.8-1
- Created the SPEC file

