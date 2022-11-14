document.querySelector("html").addEventListener("click", () => alert("clicked"));
let image: HTMLImageElement = document.querySelector("img");
image.onclick = () => alert(image.alt);