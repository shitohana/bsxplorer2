mod convert;

use crate::convert::InputFilesProps;
use bsxplorer_ci::convert::{ConvertReportType, IpcCompression};
use std::collections::HashMap;
use stylist::yew::styled_component;
use wasm_bindgen::prelude::wasm_bindgen;
use wasm_bindgen::JsValue;
use yew;
use yew::prelude::*;
use yew::props;
use yew_router::prelude::*;
#[derive(Clone, Routable, PartialEq)]
enum Route {
    #[at("/")]
    Home,
    #[at("/about")]
    About,
    #[at("/dmr")]
    Dmr,
    #[at("/convert")]
    Convert,
    #[at("/stats")]
    Stats,
    #[not_found]
    #[at("/404")]
    NotFound,
}

#[function_component(Home)]
fn home() -> Html {
    html! {
        <div>
            <h1>{ "Home Page" }</h1>
            <p>{ "Welcome to the Home page." }</p>
        </div>
    }
}

#[function_component(About)]
fn about() -> Html {
    html! {
        <div>
            <h1>{ "About Page" }</h1>
            <p>{ "Learn more About us on this page." }</p>
        </div>
    }
}

#[function_component(Stats)]
fn stats() -> Html {
    html! {
        <div>
            <h1>{ "About Page" }</h1>
            <p>{ "Learn more About us on this page." }</p>
        </div>
    }
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_name = "showDirectoryPicker")]
    fn show_directory_picker() -> JsValue;
}

#[derive(Clone, Debug, PartialEq)]
pub struct Config {
    pub input_files: Vec<String>,
    pub output_dir: String,
    pub from_type: ConvertReportType,
    pub into_type: ConvertReportType,
    pub ipc_compression: Option<IpcCompression>,
    pub threads: u32,
    pub low_memory: bool,
    pub chunk_size: u32,
    pub fasta_path: Option<String>,
    pub fai_path: Option<String>,
    pub batch_per_read: u32,
    pub batch_size: u32,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            input_files: vec![],
            output_dir: "".into(),
            from_type: ConvertReportType::Bsx,
            into_type: ConvertReportType::Bsx,
            ipc_compression: Some(IpcCompression::None),
            // In a real app you might detect the number of system threads;
            // here we default to 4 for demonstration purposes.
            threads: 4,
            low_memory: false,
            chunk_size: 20000,
            fasta_path: None,
            fai_path: None,
            batch_per_read: 8,
            batch_size: 1000000,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
enum ViewState {
    Config,
    Confirm,
    Progress,
    Done,
}

#[styled_component(Convert)]
fn convert() -> Html {
    let input_files = use_state(|| Vec::<String>::new());
    let invalid = use_state(|| HashMap::<usize, String>::new());

    let config = use_state(|| {
        props! {
            convert::Config {

            }
        }
    });

    // --- Styling with stylist ---
    let stylesheet = convert::style();

    html! {
        <div class={stylesheet}>
            {
                match config.state {
                    convert::ViewState::Config => html!{ <convert::states::ConfigView/> },
                    // convert::ViewState::Config => convert::states::config(config),
                    convert::ViewState::Confirm => convert::states::confirm(config),
                    convert::ViewState::Progress => convert::states::progress(config),
                    convert::ViewState::Done => convert::states::done(config),
            }
        }
        </div>
    }
}

#[function_component(Dmr)]
fn dmr() -> Html {
    html! {
        <div>
            <h1>{ "About Page" }</h1>
            <p>{ "Learn more About us on this page." }</p>
        </div>
    }
}

#[function_component(NavMenu)]
fn nav_menu() -> Html {
    let is_menu_open = use_state(|| false);
    let is_settings_open = use_state(|| false);

    let toggle_menu = {
        let is_menu_open = is_menu_open.clone();
        Callback::from(move |_| is_menu_open.set(!*is_menu_open))
    };

    let toggle_settings = {
        let is_settings_open = is_settings_open.clone();
        Callback::from(move |_| is_settings_open.set(!*is_settings_open))
    };

    html! {
        <nav class="navbar">
            <div class="menu-toggle" onclick={toggle_menu.clone()}>
                <div class="bar"></div>
                <div class="bar"></div>
                <div class="bar"></div>
            </div>

            // Navigation Links
            <ul class={if *is_menu_open { "nav-links open" } else { "nav-links" }}>
                <li><Link<Route> to={Route::Home}>{ "Home" }</Link<Route>></li>
                <li><Link<Route> to={Route::About}>{ "About" }</Link<Route>></li>
                <li><Link<Route> to={Route::Dmr}>{ "DMR" }</Link<Route>></li>
                <li><Link<Route> to={Route::Convert}>{ "Convert" }</Link<Route>></li>
                <li><Link<Route> to={Route::Stats}>{ "Stats" }</Link<Route>></li>
            </ul>

            // Settings Button (Top-Right Corner)
            <div class="settings-container">
                <button class="settings-btn" onclick={toggle_settings.clone()}>{ "⚙️" }</button>
                {
                    if *is_settings_open {
                        html! {
                            <div class="settings-menu">
                                <ul>
                                    <li>{ "Profile" }</li>
                                </ul>
                            </div>
                        }
                    } else {
                        html! {}
                    }
                }
            </div>
        </nav>
    }
}

#[function_component(App)]
fn app() -> Html {
    html! {
        <BrowserRouter>
            <NavMenu />
            <Switch<Route> render={move |routes| match routes {
                Route::Home => html! { <Home /> },
                Route::About => html! { <About /> },
                Route::Dmr => html! { <Dmr /> },
                Route::Convert => html! { <Convert /> },
                Route::Stats => html! { <Stats /> },
                Route::NotFound => html! { <h1>{ "404 - Page Not Found" }</h1> },
            } } />
        </BrowserRouter>
    }
}

fn main() {
    yew::Renderer::<App>::new().render();
}
