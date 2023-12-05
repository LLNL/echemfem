import numbers

from itertools import chain

from ufl.measure import Measure
from ufl.core.expr import Expr
from ufl.checks import is_true_ufl_scalar
from ufl.constantvalue import as_ufl
from ufl.domain import extract_domains


class CylindricalMeasure(Measure):
    def __init__(self, radius, *args, **kwargs):
        self.radius = radius
        super().__init__(*args, **kwargs)

    def __rmul__(self, integrand):
        """Multiply a scalar expression with measure to construct a form with
        a single integral.

        This is to implement the notation

            form = integrand * self.radius * self

        Integration properties are taken from this Measure object.

        """
        # Avoid circular imports
        from ufl.integral import Integral
        from ufl.form import Form

        # Allow python literals: 1*dx and 1.0*dx
        if isinstance(integrand, (int, float)):
            integrand = as_ufl(integrand)

        # Let other types implement multiplication with Measure if
        # they want to (to support the dolfin-adjoint TimeMeasure)
        if not isinstance(integrand, Expr):
            return NotImplemented

        # Allow only scalar integrands
        if not is_true_ufl_scalar(integrand):
            raise ValueError(
                "Can only integrate scalar expressions. The integrand is a "
                f"tensor expression with value shape {integrand.ufl_shape} and "
                f"free indices with labels {integrand.ufl_free_indices}.")

        # If we have a tuple of domain ids build the integrals one by
        # one and construct as a Form in one go.
        subdomain_id = self.subdomain_id()
        if isinstance(subdomain_id, tuple):
            return Form(list(chain(*((integrand * self.reconstruct(subdomain_id=d)).integrals()
                                     for d in subdomain_id))))

        # Check that we have an integer subdomain or a string
        # ("everywhere" or "otherwise", any more?)
        if not isinstance(subdomain_id, (str, numbers.Integral,)):
            raise ValueError("Expecting integer or string domain id.")

        # If we don't have an integration domain, try to find one in
        # integrand
        domain = self.ufl_domain()
        if domain is None:
            domains = extract_domains(integrand)
            if len(domains) == 1:
                domain, = domains
            elif len(domains) == 0:
                raise ValueError("This integral is missing an integration domain.")
            else:
                raise ValueError("Multiple domains found, making the choice of integration domain ambiguous.")

        # Otherwise create and return a one-integral form
        integral = Integral(integrand=self.radius*integrand,
                            integral_type=self.integral_type(),
                            domain=domain,
                            subdomain_id=subdomain_id,
                            metadata=self.metadata(),
                            subdomain_data=self.subdomain_data())
        return Form([integral])

    def reconstruct(self,
                    integral_type=None,
                    subdomain_id=None,
                    domain=None,
                    metadata=None,
                    subdomain_data=None):
        if subdomain_id is None:
            subdomain_id = self.subdomain_id()
        if domain is None:
            domain = self.ufl_domain()
        if metadata is None:
            metadata = self.metadata()
        if subdomain_data is None:
            subdomain_data = self.subdomain_data()
        return CylindricalMeasure(self.radius, self.integral_type(),
                                  domain=domain, subdomain_id=subdomain_id,
                                  metadata=metadata, subdomain_data=subdomain_data)
