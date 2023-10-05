using Pkg
original_active_project = Base.active_project()
try
  withenv("JULIA_PKG_SERVER" => nothing) do
    Pkg.activate("v$(VERSION.major).$(VERSION.minor)"; shared = true)
    Pkg.add("PkgAuthentication")
  end
finally
  Pkg.activate(original_active_project)
end

using PkgAuthentication
PkgAuthentication.install("juliahub.com")

Pkg.Registry.add()