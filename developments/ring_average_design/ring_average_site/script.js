const reveals = Array.from(document.querySelectorAll(".reveal"));

const io = new IntersectionObserver(
  (entries) => {
    entries.forEach((entry) => {
      if (entry.isIntersecting) {
        entry.target.classList.add("in");
        io.unobserve(entry.target);
      }
    });
  },
  { threshold: 0.12 }
);

reveals.forEach((node, idx) => {
  node.style.transitionDelay = `${Math.min(idx * 50, 220)}ms`;
  io.observe(node);
});
